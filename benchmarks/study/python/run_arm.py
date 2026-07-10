#!/usr/bin/env python3
"""benchmarks/study/python/run_arm.py -- WU-8 (leafblower-2ouc.9).

Registry-driven run-matrix driver for the Python arm (competitors.py's 8
ranked + ott-jax, and leafblower_adapter.py's 9 methods via public
leafblower.harvest()). Implements DESIGN.md Sec4/Sec5 measurement protocol:

  * ONE FRESH SUBPROCESS per (solver,problem,thread,build) cell (correct
    peak-RSS attribution) -- this script re-execs itself with `--worker`.
  * Applicability GATED by registry.json against WU-10-INSTALLED versions
    (benchmarks/study/common/run_matrix.py: build_matrix()/is_applicable()).
  * Single end-to-end `wall_time_s` per contract v2 -- the adapter already
    stops its own timer at weights-out; this driver only aggregates repeats
    (>=2 warmups discarded, median + 5/95 percentile + bootstrap CI,
    min-total-duration batching) and NEVER re-times anything itself.
  * Per-(language,package) DATA-LOADED RSS baseline subprocess (imports the
    adapter module + loads a representative problem spec, exits before
    solving) for downstream baseline subtraction (Gap I) -- reported
    alongside runs.parquet, never subtracted from the frozen raw column.
  * Hardware isolation state CAPTURED + ASSERTED (governor/turbo/pinning)
    into environment.json, warn-not-fail; run order RANDOMIZED.
  * ott-jax EXCLUDED from ranked {1,4}-thread timing (WU-OTT owns its own
    vignette) -- see run_matrix.RANKED_TIMING_EXCLUDED.
  * Pre-run frozen-tag SHA-check hook (run_matrix.assert_frozen_tag()).
  * Build variants: leafblower rows carry build in {portable,native} per
    registry.json['solvers'][...]['builds'] (WU-10-installed set only).
    This driver takes `build` as a pass-through TAG ONLY -- it never edits
    python/CMakeLists.txt or triggers a rebuild (WU-8b/leafblower-2ouc.21
    owns the portable-build swap procedure); running this driver under a
    build variant that is not currently installed will surface as adapter
    import/solve failures, not a silent mislabel.

STRICT SEPARATION (user constraint 2026-07-08): this file, and every module
it imports, reaches leafblower ONLY through leafblower_adapter.py's public
`leafblower.harvest()` wrapper. No leafblower source file is read or edited
by this WU.

Usage (repo root, single-thread BLAS forced internally per-cell):
    cd python && .venv/bin/python ../benchmarks/study/python/run_arm.py --smoke
    cd python && .venv/bin/python ../benchmarks/study/python/run_arm.py --sync-registry
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from typing import Any, Optional

_THIS_DIR = Path(__file__).resolve().parent              # benchmarks/study/python
_STUDY_DIR = _THIS_DIR.parent                               # benchmarks/study
_REPO_ROOT = _STUDY_DIR.parent.parent                         # repo root

sys.path.insert(0, str(_THIS_DIR))
sys.path.insert(0, str(_STUDY_DIR / "common"))

import run_matrix as rm  # noqa: E402

CONTRACT_KEYS = {
    "weights_ref", "iterations", "status", "converged", "error_message",
    "wall_time_s", "peak_rss_bytes",
}
STATUS_ENUM = {
    "converged", "no_conv", "infeasible", "bound_violation", "bad_arg",
    "budget", "stall", "error", "dnf",
}
RUNS_ROW_KEYS = [
    "solver", "problem", "thread", "build", "rep",
    "weights_ref", "iterations", "status", "converged", "error_message",
    "wall_time_s", "peak_rss_bytes", "trajectory_ref",
]


def _set_thread_env(thread: int) -> None:
    """MUST run before numpy/pandas/leafblower/competitors are imported in
    THIS process (CLAUDE.md single-thread-BLAS determinism rule, extended
    here to the per-cell thread sweep). Uses direct assignment (not
    setdefault) so it wins over the adapters' own `setdefault` calls."""
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[var] = str(thread)


def _self_peak_rss_bytes() -> Optional[int]:
    import resource
    try:
        return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024)
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Worker mode: ONE fresh subprocess, ONE (solver,problem,thread,build) cell.
# ---------------------------------------------------------------------------

def _resolve_adapter(solver_id: str):
    """Returns a zero/one-arg-beyond-problem callable per adapter contract:
    leafblower adapters take (problem, build); competitor adapters take
    (problem,) only (build is always 'na' for them -- registry-enforced)."""
    if solver_id.startswith("leafblower_") and solver_id.endswith("_py"):
        import leafblower_adapter
        fn = leafblower_adapter.LEAFBLOWER_PY_ADAPTERS[solver_id]
        return lambda problem, build: fn(problem, build)
    import competitors
    fn = competitors.ADAPTERS[solver_id]
    return lambda problem, build: fn(problem)


def _subsample(problem: dict, n_cap: Optional[int], seed: int) -> dict:
    """--smoke data-reduction hook: caps n via a fixed-seed random subsample
    (targets are population proportions, invariant to n, so subsampling the
    raw rows is a valid smoke-scale reduction of a real problem, not a
    different problem)."""
    if n_cap is None or len(problem["data"]) <= n_cap:
        return problem
    out = dict(problem)
    out["data"] = problem["data"].sample(n=n_cap, random_state=seed).reset_index(drop=True)
    out["design_weights"] = problem["design_weights"][
        problem["data"].sample(n=n_cap, random_state=seed).index.to_numpy()
    ] if problem["design_weights"] is not None else None
    return out


def run_worker(cell_path: str, result_path: str) -> None:
    with open(cell_path) as f:
        cell = json.load(f)

    _set_thread_env(cell["thread"])  # MUST precede numpy/pandas/adapter import below

    from problem_io import load_problem_spec  # noqa: E402 -- deferred past env-var set
    import instance_family  # noqa: E402 -- deferred past env-var set

    # Instance-family specs carry a `gen:instance_family?...` data_ref that
    # problem_io.load_problem_spec() cannot resolve on its own (raises
    # NotImplementedError). install_gen_resolver() monkey-patches the resolver;
    # it is idempotent and intercepts ONLY gen:instance_family refs, so it is a
    # no-op for the 4 static specs -- safe to install unconditionally here.
    instance_family.install_gen_resolver()

    fn = _resolve_adapter(cell["solver"])
    problem = load_problem_spec(rm.resolve_spec_path(cell["problem"]))
    problem = _subsample(problem, cell.get("n_cap"), cell.get("seed", 0))

    # Universal structural-infeasibility short-circuit: a target category with 0
    # sample observations makes the cell infeasible for EVERY solver -> record one
    # "infeasible" row, do NOT run the solver (no misleading crash/error/fake-conv).
    infeas = structural_infeasible_cats(problem)
    if infeas:
        rows = [_infeasible_row(cell, infeas)]
        summary = rm.summarize_timing([0.0], rng_seed=cell.get("seed"))
        with open(result_path, "w") as f:
            json.dump(dict(rows=rows, summary=summary), f)
        return

    warmups = int(cell["warmups"])
    min_reps = int(cell["reps"])
    min_total_duration = float(cell["min_total_duration"])
    max_reps = int(cell["max_reps"])

    timed: list[dict[str, Any]] = []
    calls = 0
    cumulative = 0.0
    while True:
        res = fn(problem, cell["build"])
        calls += 1
        if calls <= warmups:
            continue
        timed.append(res)
        cumulative += float(res["wall_time_s"])
        if len(timed) >= min_reps and cumulative >= min_total_duration:
            break
        if len(timed) >= max_reps:
            break

    weights_refs = {r["weights_ref"] for r in timed}
    if len(weights_refs) > 1:
        print(f"WARN: cell {cell['solver']}/{cell['problem']}/t{cell['thread']}/{cell['build']}: "
              f"weights_ref differs across reps ({weights_refs}) -- contract.md Sec2.1 expects a "
              f"single shared file per cell; recording as observed (not silently coerced).",
              file=sys.stderr)

    cell_peak_rss = max(
        (r["peak_rss_bytes"] for r in timed if r.get("peak_rss_bytes") is not None),
        default=_self_peak_rss_bytes(),
    )

    rows = []
    for rep_idx, r in enumerate(timed):
        assert set(r.keys()) == CONTRACT_KEYS, f"{cell['solver']}: adapter key drift {set(r.keys())}"
        assert r["status"] in STATUS_ENUM, f"{cell['solver']}: status {r['status']!r} not harmonized"
        row = dict(
            solver=cell["solver"], problem=cell["problem"], thread=cell["thread"],
            build=cell["build"], rep=rep_idx,
            weights_ref=r["weights_ref"], iterations=r["iterations"], status=r["status"],
            converged=bool(r["converged"]), error_message=r["error_message"],
            wall_time_s=float(r["wall_time_s"]), peak_rss_bytes=cell_peak_rss,
            # Per-iteration margin-error trajectory (RQ3): NOT exposed by any
            # currently-shipped adapter (WU-5/6/7's contract is the frozen
            # 7-key dict; verified via CONTRACT_KEYS assert above). This
            # column is the honest hook -- always null today, populated iff
            # a future adapter extends its return beyond the frozen 7 keys.
            trajectory_ref=None,
        )
        rows.append(row)

    summary = rm.summarize_timing([row["wall_time_s"] for row in rows], rng_seed=cell.get("seed"))
    with open(result_path, "w") as f:
        json.dump(dict(rows=rows, summary=summary), f)


def run_baseline_worker(solver_id: str, result_path: str) -> None:
    """Data-loaded RSS baseline (Gap I): imports the adapter module that
    owns `solver_id` + loads a representative problem spec, exits BEFORE
    solving. Reported alongside runs.parquet for downstream baseline
    subtraction; never mutates the raw peak_rss_bytes column."""
    _set_thread_env(1)
    t0 = time.perf_counter()
    if solver_id.startswith("leafblower_"):
        import leafblower_adapter  # noqa: F401
    else:
        import competitors  # noqa: F401
    from problem_io import load_problem_spec
    load_problem_spec(_STUDY_DIR / "spec" / "toy_inline.json")
    load_ms = time.perf_counter() - t0
    with open(result_path, "w") as f:
        json.dump(dict(solver=solver_id, load_time_s=load_ms,
                        peak_rss_bytes=_self_peak_rss_bytes()), f)


# ---------------------------------------------------------------------------
# Orchestrator mode
# ---------------------------------------------------------------------------

class _CellTimeout(RuntimeError):
    """A worker subprocess exceeded --cell-timeout (wall-clock budget). Distinct
    from a crash so orchestrate() can record a right-censored `dnf` row."""


def _spawn(args: list[str], timeout_s: float) -> None:
    import subprocess
    try:
        proc = subprocess.run([sys.executable, str(Path(__file__).resolve()), *args],
                               cwd=str(_REPO_ROOT), capture_output=True, text=True, timeout=timeout_s)
    except subprocess.TimeoutExpired as e:
        raise _CellTimeout(f"exceeded --cell-timeout={timeout_s}s: args={args}") from e
    if proc.returncode != 0:
        raise RuntimeError(f"worker subprocess failed (exit {proc.returncode}): "
                            f"args={args}\nstdout={proc.stdout}\nstderr={proc.stderr}")


def _dnf_row(cell: dict, budget_s: float) -> dict[str, Any]:
    """Right-censored did-not-finish row for a cell killed at the wall-clock
    budget -- recorded (not dropped) so a timeout is distinguishable from a cell
    that never ran. Weights are UNDEFINED for a cell that never finished, so
    weights_ref is a length-1 all-NaN sentinel (non-dangling; scoring skips
    weights for status=='dnf'). status='dnf'; wall_time_s is the budget (the
    censoring point); peak_rss unknown."""
    import numpy as np
    import pandas as pd
    weights_dir = _REPO_ROOT / "weights"
    weights_dir.mkdir(parents=True, exist_ok=True)
    path = weights_dir / f"{cell['solver']}__{cell['problem']}__t{cell['thread']}__{cell['build']}.parquet"
    pd.DataFrame({"weight": [np.nan]}).to_parquet(path)
    ref = str(path.relative_to(_REPO_ROOT))
    return dict(solver=cell["solver"], problem=cell["problem"], thread=cell["thread"],
                build=cell["build"], rep=0, weights_ref=ref, iterations=None,
                status="dnf", converged=False,
                error_message=f"exceeded --cell-timeout={budget_s}s (right-censored DNF)",
                wall_time_s=float(budget_s), peak_rss_bytes=None, trajectory_ref=None)


def structural_infeasible_cats(problem: dict) -> list:
    """Margin categories with target>0 but ZERO sample observations (no row
    carries them) -- unreachable by ANY reweighting. Extreme-skew/high-cardinality
    small-n instances realize only a subset of the declared categories. Checked on
    the MATERIALIZED (post-subsample) problem the solver would see, so the universal
    short-circuit in run_worker applies uniformly to every solver -- competitors AND
    leafblower. Mirrors run_arm.R::.structural_infeasible_cats."""
    import pandas as pd
    bad = []
    data, targets = problem["data"], problem["targets"]
    for m in problem["margins"]:
        present = set(map(str, pd.unique(data[m])))
        bad += [f"{m}.{k}" for k, v in targets[m].items() if v > 0 and str(k) not in present]
    return bad


def _infeasible_row(cell: dict, cats: list) -> dict[str, Any]:
    """Row synthesized when a cell is structurally infeasible: no solver runs;
    weights undefined -> length-1 all-NaN sentinel (scoring skips), status=
    'infeasible'. Mirrors _dnf_row's sentinel; wall_time_s=0 (nothing ran)."""
    import numpy as np
    import pandas as pd
    weights_dir = _REPO_ROOT / "weights"
    weights_dir.mkdir(parents=True, exist_ok=True)
    path = weights_dir / f"{cell['solver']}__{cell['problem']}__t{cell['thread']}__{cell['build']}.parquet"
    pd.DataFrame({"weight": [np.nan]}).to_parquet(path)
    ref = str(path.relative_to(_REPO_ROOT))
    more = f" (+{len(cats) - 3} more)" if len(cats) > 3 else ""
    return dict(solver=cell["solver"], problem=cell["problem"], thread=cell["thread"],
                build=cell["build"], rep=0, weights_ref=ref, iterations=None,
                status="infeasible", converged=False,
                error_message=f"structurally infeasible: target>0 but 0 observations for "
                              f"{', '.join(cats[:3])}{more}",
                wall_time_s=0.0, peak_rss_bytes=None, trajectory_ref=None)


def _crash_row(cell: dict, reason: str) -> dict[str, Any]:
    """Worker died without a result (OOM-kill / hard crash, non-timeout): recorded
    as status='error' (NOT dropped to cell_failures) so a crashed cell does not
    VANISH -- the WU-9 no-selective-reporting guarantee, extending _dnf_row's
    treatment from timeouts to crashes. Weights undefined -> NaN sentinel; timing
    unknown -> NaN. Mirrors run_arm.R::.crash_row."""
    import numpy as np
    import pandas as pd
    weights_dir = _REPO_ROOT / "weights"
    weights_dir.mkdir(parents=True, exist_ok=True)
    path = weights_dir / f"{cell['solver']}__{cell['problem']}__t{cell['thread']}__{cell['build']}.parquet"
    pd.DataFrame({"weight": [np.nan]}).to_parquet(path)
    return dict(solver=cell["solver"], problem=cell["problem"], thread=cell["thread"],
                build=cell["build"], rep=0, weights_ref=str(path.relative_to(_REPO_ROOT)),
                iterations=None, status="error", converged=False,
                error_message=f"worker died without result: {reason[:160]}",
                wall_time_s=float("nan"), peak_rss_bytes=None, trajectory_ref=None)


def _available_adapters() -> set:
    """Solver ids that resolve to a runnable Python adapter this arm: the wrapped
    competitors plus the 9 leafblower_*_py methods. A registry solver absent here
    (e.g. `samplics`, installed but deliberately unwrapped -- svy is its successor)
    must be gated out of the matrix, else _resolve_adapter KeyErrors on it."""
    import competitors
    import leafblower_adapter
    return set(competitors.ADAPTERS) | set(leafblower_adapter.LEAFBLOWER_PY_ADAPTERS)


def orchestrate(opts: argparse.Namespace) -> dict[str, Any]:
    registry = rm.load_registry()
    # Full matrix problem set: the 4 static specs UNION the WU-3 instance family
    # (32 instances). The instance family is the only arena declaring
    # ot/newton_kl objective_families, so it is what gives the OT / Newton-KL
    # solvers their applicable problems. This also materializes the generated
    # spec/instance_family/<id>.json files the worker later loads.
    # regenerate=False: use the COMMITTED (frozen) instance specs -- regenerating
    # would overwrite them and dirty the frozen runnable tree (WU-9). The committed
    # spec/instance_family/*.json ARE the frozen truth; regen only on a missing dir.
    problem_specs = rm.load_all_problem_specs(regenerate=False)
    hyperparams = rm.load_hyperparams()
    available = _available_adapters()

    if opts.sync_registry:
        registry = rm.compute_applicable_problems(registry, problem_specs, available=available, arm="python")
        with open(_STUDY_DIR / "registry.json", "w") as f:
            json.dump(registry, f, indent=2)
            f.write("\n")
        # --sync-registry is SYNC-ONLY: it refreshes applicable_problems and
        # exits, never launching a run. A bounded --smoke run may be combined
        # with it (sync then smoke); a bare --sync-registry must NOT fall
        # through to the full production matrix (that is WU-11, post-freeze).
        if not opts.smoke:
            return dict(synced=True, n_solvers=len(registry["solvers"]),
                        registry_path=str(_STUDY_DIR / "registry.json"))

    tag_status = rm.assert_frozen_tag()
    # WU-9: a scored/timed run must execute against the frozen runnable tree. When
    # --assert-runnable-tag is passed (WU-11 scored launcher), hard-stop on a
    # dirty/drifted benchmarks/study tree. The rehearsal omits it (unfrozen by design).
    if opts.assert_runnable_tag:
        rm.assert_runnable_frozen(opts.assert_runnable_tag)
        print(f"WU-9 runnable-tree gate OK: tree matches signed tag {opts.assert_runnable_tag!r}",
              file=sys.stderr)
    env_info = rm.capture_environment(pin_core=opts.pin_core)

    # --smoke bounds the instance family to instances with n <= opts.n: a huge-n
    # synthetic instance would burn minutes in the pure-Python generator only to
    # be subsampled down to opts.n rows. Static specs (instance_family_n is None)
    # are always kept -- they load from file/pkg and are subsampled cheaply.
    matrix_specs = problem_specs
    if opts.smoke:
        matrix_specs = {
            pid: spec for pid, spec in problem_specs.items()
            if (rm.instance_family_n(spec) is None or rm.instance_family_n(spec) <= opts.n)
        }

    threads = tuple(int(t) for t in opts.threads.split(","))
    cells = rm.build_matrix(registry, matrix_specs, arm="python", threads=threads, available=available,
                             rng_seed=opts.seed)
    if opts.solvers:
        import re
        pat = re.compile(opts.solvers)
        cells = [c for c in cells if pat.search(c["solver"])]
    if opts.problems:
        import re
        ppat = re.compile(opts.problems)
        cells = [c for c in cells if ppat.search(c["problem"])]
    if opts.max_cells:
        cells = cells[: opts.max_cells]

    out_dir = Path(opts.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    work_dir = out_dir / "_cells"
    work_dir.mkdir(parents=True, exist_ok=True)

    warmups = 0 if opts.smoke else opts.warmups
    reps = opts.reps
    min_total_duration = 0.0 if opts.smoke else opts.min_total_duration
    n_cap = opts.n if opts.smoke else None

    all_rows: list[dict[str, Any]] = []
    cell_failures: list[dict[str, Any]] = []
    for i, cell in enumerate(cells):
        cell_full = dict(cell, warmups=warmups, reps=reps, max_reps=opts.max_reps,
                          min_total_duration=min_total_duration, n_cap=n_cap, seed=opts.seed)
        cell_path = work_dir / f"cell_{i}.json"
        result_path = work_dir / f"result_{i}.json"
        with open(cell_path, "w") as f:
            json.dump(cell_full, f)
        try:
            _spawn(["--worker", "--cell", str(cell_path), "--result-out", str(result_path)],
                   timeout_s=opts.cell_timeout)
            with open(result_path) as f:
                out = json.load(f)
            all_rows.extend(out["rows"])
        except _CellTimeout:
            # Right-censored DNF: record a real row, do NOT drop (else a timeout
            # is indistinguishable from 'never ran').
            all_rows.append(_dnf_row(cell, opts.cell_timeout))
            print(f"DNF: cell {cell['solver']}/{cell['problem']}/t{cell['thread']}/{cell['build']} "
                  f"exceeded {opts.cell_timeout}s budget (right-censored)", file=sys.stderr)
        except Exception as e:  # noqa: BLE001 -- one cell's crash must not abort the matrix
            # Worker died without a result (OOM-kill / hard crash). Record a real
            # error row -- do NOT drop, else a crashed cell VANISHES (WU-9 forbids
            # selective reporting).
            all_rows.append(_crash_row(cell, str(e)))
            cell_failures.append(dict(cell=cell, error=str(e)))
            print(f"CRASH: cell {cell['solver']}/{cell['problem']} died -> recorded as error row: {e}",
                  file=sys.stderr)

    baseline_reps: dict[str, str] = {}
    for cell in cells:
        pkg = rm.solver_package(cell["solver"], hyperparams)
        baseline_reps.setdefault(pkg, cell["solver"])
    baselines = []
    for pkg, solver_id in baseline_reps.items():
        result_path = work_dir / f"baseline_{pkg}.json"
        try:
            _spawn(["--baseline-rss", "--solver", solver_id, "--result-out", str(result_path)],
                   timeout_s=opts.cell_timeout)
            with open(result_path) as f:
                b = json.load(f)
            b["package"] = pkg
            baselines.append(b)
        except Exception as e:  # noqa: BLE001
            print(f"WARN: baseline for package {pkg!r} (via {solver_id}) failed: {e}", file=sys.stderr)

    import pandas as pd
    runs_df = pd.DataFrame(all_rows, columns=RUNS_ROW_KEYS) if all_rows else pd.DataFrame(columns=RUNS_ROW_KEYS)
    if not runs_df.empty:
        runs_df["thread"] = runs_df["thread"].astype("int32")
        runs_df["build"] = runs_df["build"].astype("category")
        runs_df["rep"] = runs_df["rep"].astype("int32")
        runs_df["status"] = runs_df["status"].astype(
            pd.CategoricalDtype(categories=sorted(STATUS_ENUM)))
    runs_path = out_dir / "runs.parquet"
    runs_df.to_parquet(runs_path, index=False)

    with open(out_dir / "environment.json", "w") as f:
        json.dump(dict(env_info, frozen_tag=tag_status, baselines=baselines,
                        n_cells=len(cells), n_cell_failures=len(cell_failures)), f, indent=2)

    return dict(n_cells=len(cells), n_rows=len(all_rows), n_failures=len(cell_failures),
                runs_path=str(runs_path), cell_failures=cell_failures)


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--worker", action="store_true", help="internal: run one cell in a fresh subprocess")
    p.add_argument("--baseline-rss", action="store_true", help="internal: data-loaded RSS baseline subprocess")
    p.add_argument("--cell", type=str, help="worker mode: path to cell JSON")
    p.add_argument("--solver", type=str, help="baseline-rss mode: solver id to import")
    p.add_argument("--result-out", type=str, help="worker/baseline-rss mode: path to write result JSON")
    p.add_argument("--smoke", action="store_true", help="n=5000 subsample, reps=2, no warmups, no duration floor")
    p.add_argument("--n", type=int, default=5000, help="--smoke row-count cap")
    p.add_argument("--reps", type=int, default=2, help="minimum timed reps per cell")
    p.add_argument("--warmups", type=int, default=2, help="warmup calls discarded per cell (ignored under --smoke)")
    p.add_argument("--max-reps", type=int, default=200, help="safety cap on timed reps (min-duration batching)")
    p.add_argument("--min-total-duration", type=float, default=0.5,
                   help="min cumulative timed seconds per cell before stopping (batching floor)")
    p.add_argument("--threads", type=str, default="1,4", help="comma-separated thread sweep")
    p.add_argument("--pin-core", type=int, default=None, help="core pinned via taskset/numactl by the caller")
    p.add_argument("--seed", type=int, default=None, help="RNG seed for cell-order shuffle + bootstrap")
    p.add_argument("--out", type=str, default=str(_STUDY_DIR / "results"), help="output directory")
    p.add_argument("--sync-registry", action="store_true",
                   help="recompute + write back registry.json's applicable_problems before running")
    p.add_argument("--solvers", type=str, default=None, help="regex filter on solver id")
    p.add_argument("--problems", type=str, default=None, help="regex filter on problem id (shard by n-tier)")
    p.add_argument("--assert-runnable-tag", type=str, default=None,
                   help="WU-9: hard-stop unless benchmarks/study is clean + matches this signed runnable-tree tag (scored runs only)")
    p.add_argument("--max-cells", type=int, default=None, help="cap number of cells (debugging/CI)")
    p.add_argument("--cell-timeout", type=float, default=120.0,
                   help="max seconds per fresh-subprocess cell before it is killed and logged as a "
                        "cell_failure (a hung competitor must not stall the whole matrix)")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_argparser().parse_args(argv)
    if args.worker:
        run_worker(args.cell, args.result_out)
        return 0
    if args.baseline_rss:
        run_baseline_worker(args.solver, args.result_out)
        return 0
    result = orchestrate(args)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
