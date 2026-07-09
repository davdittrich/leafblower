"""benchmarks/study/python/test_adapter_ott.py -- WU-OTT (leafblower-2ouc.18).

Tests for the dedicated ott-jax (JAX/XLA) adapter carved out of
competitors.py. Strict separation: this file imports nothing from
leafblower's own src/, r_bridge.cpp, R/, python/leafblower/.

Run with (single-thread BLAS, venv python, per CLAUDE.md determinism rule):
  cd python && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    .venv/bin/python -m pytest ../benchmarks/study/python/test_adapter_ott.py -q
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import sys
from pathlib import Path

import numpy as np

_THIS_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_THIS_DIR))

import adapter_ott  # noqa: E402
import competitors  # noqa: E402
import test_competitors as _tc  # noqa: E402  (reuse the SAME hand-solved IPF golden as pot_sinkhorn)


# ---------------------------------------------------------------------------
# 1. x64 actually in effect (the load-bearing DoD: a test that would FAIL
#    under JAX's float32 default -- proves jax_enable_x64 genuinely applied,
#    not merely requested).
# ---------------------------------------------------------------------------

def test_x64_actually_in_effect():
    import jax.numpy as jnp

    z = jnp.zeros(1, dtype=jnp.float64)
    assert z.dtype == jnp.float64, (
        f"jax_enable_x64 not in effect: dtype={z.dtype} (would be float32 "
        "under JAX's default -- adapter_ott.py's module-top "
        "jax.config.update('jax_enable_x64', True) failed to apply)"
    )
    # Numeric proof beyond the dtype label: delta is representable
    # distinctly at float64 precision (eps ~2.22e-16) but rounds away at
    # float32 precision (eps ~1.19e-7). If JAX silently downcast storage
    # to float32 despite the float64 dtype label, val64 would collapse to
    # numpy's float32-rounded reference instead of its float64 one.
    delta = 1e-9
    val64 = float(jnp.asarray(1.0 + delta, dtype=jnp.float64))
    ref32 = float(np.float32(1.0 + delta))
    ref64 = float(np.float64(1.0 + delta))
    assert val64 == ref64, "float64 precision lost -- silently running in float32"
    assert val64 != ref32, "delta collapsed to the float32-rounded value"


def test_jax_enable_x64_config_flag_is_true():
    import jax

    assert jax.config.jax_enable_x64 is True


# ---------------------------------------------------------------------------
# 2. XLA_FLAGS recorded / disclosed (DESIGN.md Section 2 disclosure
#    requirement).
# ---------------------------------------------------------------------------

def test_xla_config_snapshot_discloses_flags_and_x64():
    snap = adapter_ott.xla_config_snapshot()
    assert snap["xla_flags"] == os.environ.get("XLA_FLAGS", "")
    assert "xla_cpu_multi_thread_eigen=false" in snap["xla_flags"]
    assert snap["jax_enable_x64"] is True
    assert isinstance(snap["jax_local_device_count"], int) and snap["jax_local_device_count"] >= 1
    assert isinstance(snap["jax_default_backend"], str) and snap["jax_default_backend"]
    # This WU's verified finding: env-var-only XLA thread pinning is NOT
    # effective in jaxlib 0.10.2 -- must be disclosed honestly, not silently
    # claimed as working (module docstring; DESIGN.md Section 2 rationale
    # for excluding ott-jax from ranked timing).
    assert snap["intra_op_threads_env_pinning_verified_effective"] is False
    assert "intra_op_parallelism_threads" in snap["note"]


def test_xla_flags_set_before_jax_import():
    # XLA_FLAGS must be a real process env var (not merely returned by the
    # snapshot function) since XLA reads it at process/backend init time,
    # before any Python-level config call could apply it.
    assert "XLA_FLAGS" in os.environ
    assert "xla_cpu_multi_thread_eigen=false" in os.environ["XLA_FLAGS"]


# ---------------------------------------------------------------------------
# 3. coupling -> weight validated vs the 2-margin IPF golden (SAME golden
#    basis as pot_sinkhorn/pot_greenkhorn, per WU-OTT spec -- imported
#    directly from test_competitors.py so both OT arms pin to the
#    identical IPF reference, not a re-derivation that could drift).
# ---------------------------------------------------------------------------

def test_ott_jax_matches_ipf_golden_same_basis_as_pot_sinkhorn():
    problem = _tc._unbounded_problem("ott-adapter-ipf-golden")
    res = adapter_ott.run_ott_jax_sinkhorn(problem)
    _tc._assert_contract_shape(res, len(problem["data"]))
    assert res["status"] == "converged", res["error_message"]
    assert res["converged"] is True

    expected = _tc._expected_ot_weights(problem)
    got = _tc._read_weights(res)
    np.testing.assert_allclose(got, expected, rtol=1e-3, atol=1e-6)

    # Cross-check: the carved adapter must agree with pot_sinkhorn on the
    # identical golden problem (both reduce to the same IPF fixed point).
    pot_res = competitors.run_pot_sinkhorn(_tc._unbounded_problem("pot-sinkhorn-cross-check"))
    pot_w = _tc._read_weights(pot_res)
    np.testing.assert_allclose(got, pot_w, rtol=1e-3, atol=1e-6)


# ---------------------------------------------------------------------------
# 4. K != 2 -> bad_arg
# ---------------------------------------------------------------------------

def test_k_not_two_is_bad_arg():
    problem = _tc._single_margin_problem("ott-adapter-bad-arg-k")
    res = adapter_ott.run_ott_jax_sinkhorn(problem)
    _tc._assert_contract_shape(res, len(problem["data"]))
    assert res["status"] == "bad_arg"
    assert res["converged"] is False
    w = _tc._read_weights(res)
    assert np.all(np.isnan(w))


# ---------------------------------------------------------------------------
# 5. competitors.py re-export still resolves (wiring constraint: existing
#    consumers -- run_arm.py's ADAPTERS[solver_id] lookup, test_competitors.py
#    -- must keep working after the carve-out).
# ---------------------------------------------------------------------------

def test_competitors_reexport_resolves_and_matches_adapter():
    assert competitors.ADAPTERS["ott_jax_sinkhorn"] is competitors.run_ott_jax_sinkhorn
    assert competitors.run_ott_jax_sinkhorn.solver_id == "ott_jax_sinkhorn"

    problem = _tc._unbounded_problem("ott-reexport-golden")
    res = competitors.run_ott_jax_sinkhorn(problem)
    _tc._assert_sane_converged(res, problem)


# ---------------------------------------------------------------------------
# 6. Scaling vignette runs + emits disclosed config. Standalone diagnostic,
#    excluded from ranked {1,4}-thread timing (WU-8 run matrix).
# ---------------------------------------------------------------------------

def test_ott_scaling_vignette_runs_and_discloses_config():
    out = adapter_ott.ott_scaling_vignette([20, 50])
    assert out["excluded_from_ranked_timing"] is True
    assert out["config"]["jax_enable_x64"] is True
    assert "xla_cpu_multi_thread_eigen=false" in out["config"]["xla_flags"]

    assert [row["n"] for row in out["rows"]] == [20, 50]
    for row in out["rows"]:
        assert row["status"] == "converged", row
        assert row["converged"] is True
        assert isinstance(row["wall_time_s"], float) and row["wall_time_s"] > 0.0


def test_ott_scaling_vignette_pin_single_core_restores_affinity():
    if not hasattr(os, "sched_getaffinity"):
        return  # platform without POSIX affinity API -- nothing to verify
    before = os.sched_getaffinity(0)
    out = adapter_ott.ott_scaling_vignette([20], pin_single_core=True)
    after = os.sched_getaffinity(0)
    assert after == before, "ott_scaling_vignette must restore the prior CPU affinity mask"
    assert out["config"]["pinned_single_core"] is True
    assert out["config"]["pinned_affinity_mask"] and len(out["config"]["pinned_affinity_mask"]) == 1
