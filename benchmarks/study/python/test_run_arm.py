"""benchmarks/study/python/test_run_arm.py -- WU-8 (leafblower-2ouc.9).

TDD tests for run_arm.py (orchestrator + worker driver) and its shared
common/run_matrix.py helpers. Covers: applicability gating against the real
registry.json + real problem specs, cell-matrix enumeration (thread sweep,
build variants, ott-jax exclusion, randomization), the frozen-tag hook
against the real repo state, hardware-isolation capture shape, adapter
resolution (including the disclosed cvxr_reference/samplics registry
entries that have no adapter implementation), and two REAL end-to-end
integration tests: a single worker cell (toy_inline + leafblower_oris_py)
and a bounded --smoke run, both schema-checked against runs_schema.json.

Single-thread BLAS is forced before numpy import (CLAUDE.md rule).
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import pandas as pd  # noqa: E402
import pytest  # noqa: E402

_THIS_DIR = Path(__file__).resolve().parent
_STUDY_DIR = _THIS_DIR.parent
_REPO_ROOT = _STUDY_DIR.parent.parent

sys.path.insert(0, str(_THIS_DIR))
sys.path.insert(0, str(_STUDY_DIR / "common"))

import run_matrix as rm  # noqa: E402
import run_arm  # noqa: E402


@pytest.fixture(scope="module")
def registry():
    # registry.json on disk ships with applicable_problems == [] for every
    # entry (WU-C did not populate it -- it is explicitly outside the
    # frozen perimeter and this WU is authorized to recompute it). Tests
    # that need cell enumeration operate on the computed copy; tests that
    # check the raw pin fields (arm/families/bounds/K_max/builds) can use
    # either since compute_applicable_problems() never touches those.
    raw = rm.load_registry()
    return rm.compute_applicable_problems(raw, rm.load_problem_specs())


@pytest.fixture(scope="module")
def problem_specs():
    return rm.load_problem_specs()


@pytest.fixture(scope="module")
def hyperparams():
    return rm.load_hyperparams()


# ---------------------------------------------------------------------------
# Applicability gating (registry.json vs real problem specs)
# ---------------------------------------------------------------------------

def test_ipfn_unbounded_only_gates_out_on_bounded_stepstone(registry, problem_specs):
    # ipfn: families=['kl'], bounds='unbounded' (frozen pin). stepstone_bounded
    # has objective_families ['kl','chi2','logit'] but bounds.max=5 (bounded)
    # -- family matches but bounds incompatibility must gate it out.
    entry = registry["solvers"]["ipfn"]
    assert entry["bounds"] == "unbounded"
    assert not rm.is_applicable(entry, problem_specs["stepstone_bounded"])
    assert rm.is_applicable(entry, problem_specs["stepstone_unbounded"])


def test_python_arm_gates_out_pkg_data_ref_problems(registry, problem_specs):
    # canonical_survey_apistrat has data_ref='pkg:survey::apistrat', an
    # R-package dataset. problem_io.py's _PKG_LOADERS is empty (no Python
    # loader registered) -- scheduling this spec on any python-arm solver
    # is guaranteed to fail with NotImplementedError at solve time, so
    # is_applicable() must gate it out regardless of family/bounds/K match.
    # Regression test for a real full-scale --smoke run failure: weightipy
    # (python, families=['kl'], bounds='both') otherwise passes every other
    # gate against canonical_survey_apistrat (kl, bounds.max=None, K=1).
    spec = problem_specs["canonical_survey_apistrat"]
    assert spec["data_ref"] == "pkg:survey::apistrat"
    weightipy = registry["solvers"]["weightipy"]
    assert weightipy["arm"] == "python"
    assert not rm.is_applicable(weightipy, spec)
    # R arm has no such restriction (problem_io.R resolves pkg: generically
    # via requireNamespace()) -- a family/bounds/K-compatible R solver must
    # stay applicable. ebal: arm='R', families=['kl'], bounds='unbounded' --
    # same family/bounds shape as weightipy above, only the arm differs.
    ebal = registry["solvers"]["ebal"]
    assert ebal["arm"] == "R"
    assert rm.is_applicable(ebal, spec)


def test_leafblower_chebyshev_minimax_family_applicable_to_stepstone_only(registry, problem_specs):
    # 4-family taxonomy realignment (2026-07-09): leafblower_chebyshev_{r,py}
    # declare families=['minimax'], bounds='both'. stepstone_bounded/
    # stepstone_unbounded now declare 'minimax' in objective_families, so
    # chebyshev is applicable to both; toy_inline/canonical_survey_apistrat
    # stay kl-only and still gate it out.
    entry = registry["solvers"]["leafblower_chebyshev_py"]
    assert rm.is_applicable(entry, problem_specs["stepstone_bounded"])
    assert rm.is_applicable(entry, problem_specs["stepstone_unbounded"])
    assert not rm.is_applicable(entry, problem_specs["toy_inline"])
    assert not rm.is_applicable(entry, problem_specs["canonical_survey_apistrat"])


def test_ott_jax_K_max_gates_out_high_K_stepstone(registry, problem_specs):
    # ott_jax_sinkhorn is a KL-family algorithm (method_class='sinkhorn');
    # 'ot' is retired as an objective_families value (4-family taxonomy
    # realignment, 2026-07-09) -- the registry pin and synthetic probe specs
    # below use 'kl' to isolate the K_max gating branch from family matching.
    entry = registry["solvers"]["ott_jax_sinkhorn"]
    assert entry["families"] == ["kl"]
    assert entry["K_max"] == 2
    assert problem_specs["stepstone_bounded"]["K"] == 9
    # stepstone_bounded now shares the 'kl' family with ott_jax_sinkhorn --
    # K_max=2 vs K=9 alone gates it out (not a family mismatch anymore).
    synthetic_kl_low_k = dict(objective_families=["kl"], bounds={"min": 0, "max": None}, K=2)
    synthetic_kl_high_k = dict(objective_families=["kl"], bounds={"min": 0, "max": None}, K=9)
    assert rm.is_applicable(entry, synthetic_kl_low_k)
    assert not rm.is_applicable(entry, synthetic_kl_high_k)
    assert not rm.is_applicable(entry, problem_specs["stepstone_bounded"])


def test_available_specs_declare_exactly_four_family_taxonomy(problem_specs):
    # 4-family taxonomy realignment (2026-07-09): available specs now declare
    # exactly the authoritative {kl,chi2,logit,minimax} union -- stepstone_
    # bounded/unbounded added 'minimax' (Chebyshev L-infinity family); 'ot'
    # and 'newton_kl' are retired as objective_families values (they are
    # method_class/algorithm-class values in registry.json now, distinct
    # from the objective family axis).
    seen_families = set()
    for spec in problem_specs.values():
        seen_families |= set(spec["objective_families"])
    assert seen_families == {"kl", "chi2", "logit", "minimax"}
    assert "ot" not in seen_families
    assert "newton_kl" not in seen_families


def test_compute_applicable_problems_preserves_frozen_pin_fields(registry, problem_specs):
    out = rm.compute_applicable_problems(registry, problem_specs)
    for sid, entry in registry["solvers"].items():
        out_entry = out["solvers"][sid]
        for k in ("arm", "families", "bounds", "K_max", "builds"):
            assert out_entry[k] == entry[k], f"{sid}.{k} frozen pin field mutated"
        assert isinstance(out_entry["applicable_problems"], list)


# ---------------------------------------------------------------------------
# Cell-matrix enumeration
# ---------------------------------------------------------------------------

def test_build_matrix_python_arm_excludes_ott_jax(registry, problem_specs):
    cells = rm.build_matrix(registry, problem_specs, arm="python", rng_seed=1)
    assert all(c["solver"] != "ott_jax_sinkhorn" for c in cells)


def test_build_matrix_covers_thread_sweep_and_only_python_arm(registry, problem_specs):
    cells = rm.build_matrix(registry, problem_specs, arm="python", threads=(1, 4), rng_seed=1)
    assert cells, "expected at least one applicable python-arm cell"
    threads_seen = {c["thread"] for c in cells}
    assert threads_seen == {1, 4}
    solvers_seen = {c["solver"] for c in cells}
    r_solvers = {sid for sid, e in registry["solvers"].items() if e["arm"] == "R"}
    assert solvers_seen.isdisjoint(r_solvers)


def test_build_matrix_leafblower_rows_carry_own_build_tag(registry, problem_specs):
    cells = rm.build_matrix(registry, problem_specs, arm="python", rng_seed=1)
    lbw_cells = [c for c in cells if c["solver"].startswith("leafblower_")]
    assert lbw_cells and all(c["build"] == "portable" for c in lbw_cells)


def test_build_matrix_randomizes_order(registry, problem_specs):
    cells_a = rm.build_matrix(registry, problem_specs, arm="python", rng_seed=1)
    cells_b = rm.build_matrix(registry, problem_specs, arm="python", rng_seed=2)
    assert len(cells_a) == len(cells_b) > 1
    assert cells_a != cells_b  # different seeds -- different orders (astronomically unlikely to tie)


# ---------------------------------------------------------------------------
# Frozen-tag hook + hardware-isolation capture
# ---------------------------------------------------------------------------

def test_assert_frozen_tag_no_drift_against_real_repo():
    status = rm.assert_frozen_tag()  # raises RuntimeError internally on drift
    assert status["ok"] is True


def test_capture_environment_shape():
    env = rm.capture_environment()
    state = env["hardware_state_for_timed_runs"]
    for k in ("cpu_governor_all_cores", "cpu_governor_conformant",
              "turbo_boost_enabled", "turbo_boost_conformant", "core_pinning_active"):
        assert k in state
    assert isinstance(env["isolated"], bool)


def test_summarize_timing_single_rep_and_multi_rep():
    single = rm.summarize_timing([0.1])
    assert single["n_reps"] == 1
    multi = rm.summarize_timing([0.10, 0.12, 0.11, 0.13, 0.09], rng_seed=7)
    assert multi["p05"] <= multi["median"] <= multi["p95"]
    assert multi["boot_ci_lo"] <= multi["boot_ci_hi"]


# ---------------------------------------------------------------------------
# Adapter resolution (including disclosed registry/adapter gaps)
# ---------------------------------------------------------------------------

def test_resolve_adapter_leafblower_python():
    fn = run_arm._resolve_adapter("leafblower_oris_py")
    assert callable(fn)


def test_resolve_adapter_competitor():
    fn = run_arm._resolve_adapter("ipfn")
    assert callable(fn)


def test_resolve_adapter_samplics_has_no_implementation():
    # registry.json lists a 'samplics' python competitor entry, but
    # competitors.py's ADAPTERS dict (frozen) has no 'samplics' key --
    # disclosed gap: this WU does not add adapters (out of file scope).
    import competitors
    assert "samplics" not in competitors.ADAPTERS
    with pytest.raises(KeyError):
        run_arm._resolve_adapter("samplics")


# ---------------------------------------------------------------------------
# Real end-to-end integration: one worker cell
# ---------------------------------------------------------------------------

def test_worker_end_to_end_toy_inline_oris(tmp_path):
    cell = dict(solver="leafblower_oris_py", problem="toy_inline", thread=1, build="portable",
                warmups=1, reps=2, max_reps=5, min_total_duration=0.0, n_cap=None, seed=0)
    cell_path = tmp_path / "cell.json"
    result_path = tmp_path / "result.json"
    cell_path.write_text(json.dumps(cell))

    proc = subprocess.run(
        [sys.executable, str(_THIS_DIR / "run_arm.py"), "--worker",
         "--cell", str(cell_path), "--result-out", str(result_path)],
        cwd=str(_REPO_ROOT), capture_output=True, text=True,
    )
    assert proc.returncode == 0, f"worker failed: stdout={proc.stdout}\nstderr={proc.stderr}"
    out = json.loads(result_path.read_text())
    rows = out["rows"]
    assert len(rows) == 2  # reps=2, warmups discarded
    for i, row in enumerate(rows):
        assert set(row.keys()) == set(run_arm.RUNS_ROW_KEYS)
        assert row["solver"] == "leafblower_oris_py"
        assert row["problem"] == "toy_inline"
        assert row["thread"] == 1
        assert row["build"] == "portable"
        assert row["rep"] == i
        assert row["status"] in run_arm.STATUS_ENUM
        assert row["wall_time_s"] > 0
        assert (_REPO_ROOT / row["weights_ref"]).exists()
    # weights_ref identical across reps (deterministic solve, contract.md Sec2.1).
    # NOTE: the frozen leafblower_adapter.py names the weights file after the
    # bare *method* ("leafblower_oris"), not the registry solver id
    # ("leafblower_oris_py") -- a frozen-adapter naming convention, not this
    # driver's choice; see leafblower_adapter.py's run_leafblower() solver_id.
    assert len({r["weights_ref"] for r in rows}) == 1
    assert Path(rows[0]["weights_ref"]).name == "leafblower_oris__toy_inline__t1__portable.parquet"


def test_baseline_rss_worker_end_to_end(tmp_path):
    result_path = tmp_path / "baseline.json"
    proc = subprocess.run(
        [sys.executable, str(_THIS_DIR / "run_arm.py"), "--baseline-rss",
         "--solver", "leafblower_oris_py", "--result-out", str(result_path)],
        cwd=str(_REPO_ROOT), capture_output=True, text=True,
    )
    assert proc.returncode == 0, f"baseline worker failed: stdout={proc.stdout}\nstderr={proc.stderr}"
    out = json.loads(result_path.read_text())
    assert out["solver"] == "leafblower_oris_py"
    assert out["peak_rss_bytes"] > 0
    assert out["load_time_s"] >= 0


# ---------------------------------------------------------------------------
# Real end-to-end integration: bounded --smoke run, schema-checked
# ---------------------------------------------------------------------------

def test_smoke_run_produces_schema_conformant_parquet(tmp_path):
    out_dir = tmp_path / "results"
    proc = subprocess.run(
        [sys.executable, str(_THIS_DIR / "run_arm.py"), "--smoke", "--sync-registry",
         "--solvers", "^(leafblower_oris_py|ipfn)$", "--threads", "1",
         "--out", str(out_dir), "--seed", "3"],
        cwd=str(_REPO_ROOT), capture_output=True, text=True,
    )
    assert proc.returncode == 0, f"smoke run failed: stdout={proc.stdout}\nstderr={proc.stderr}"

    runs_path = out_dir / "runs.parquet"
    assert runs_path.exists()
    df = pd.read_parquet(runs_path)
    assert len(df) > 0, "smoke run produced zero rows -- registry gating or adapter resolution broke"
    with open(_STUDY_DIR / "spec" / "runs_schema.json") as f:
        schema = json.load(f)
    required = schema["items"]["required"]
    for col in required:
        assert col in df.columns, f"runs.parquet missing required column {col}"
    assert set(df["status"].dropna().unique()) <= run_arm.STATUS_ENUM
    assert set(df["thread"].unique()) <= {1, 4}
    assert set(df["build"].unique()) <= {"portable", "native", "na"}
    assert (df["wall_time_s"] > 0).all()
    assert df["rep"].min() == 0

    env_path = out_dir / "environment.json"
    assert env_path.exists()
    env = json.loads(env_path.read_text())
    assert "hardware_state_for_timed_runs" in env
    assert "frozen_tag" in env
    assert "baselines" in env
