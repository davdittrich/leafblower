"""benchmarks/study/python/leafblower_adapter.py -- WU-7 (leafblower-2ouc.8).

Uniform-contract Python adapters wrapping leafblower's OWN 9 solver methods
(oris, oris_soft, raking, sinkhorn, greenkhorn, chebyshev, greg, logit,
newton_kl), timed under TWO build variants (portable / native).

STRICT SEPARATION (user constraint, 2026-07-08): this file calls leafblower
ONLY through its PUBLIC `leafblower.harvest()` entry point. No leafblower
source file (src/, r_bridge.cpp, python/leafblower/, pyproject.toml) is read
or edited here, and there is no solve-only entry point -- leafblower is
timed end-to-end (groupby/cell-build inside the timer) exactly like a
whole-unit competitor (contract.md Section 1.1).

Contract v2 (benchmarks/study/spec/contract.md): every `run_leafblower()`
call below returns EXACTLY
    {weights_ref, iterations, status, converged, error_message,
     wall_time_s, peak_rss_bytes}
`status` is leafblower's OWN native RK_* status code mapped to the 8-value
harmonized enum (status_enum.json), with a harness-side bound_violation
reclassification layered on top. `converged` is ALWAYS independently
recomputed here from common/metrics.py's margin_stats() against
problem["tol"] -- never the solver's self-reported status (contract.md
Section 2.4). Python's harvest() cleanly distinguishes infeasible
(RuntimeError, status==2) from bad_arg (ValueError, status==3 and all
pre-validation errors) by exception TYPE.

`build` (portable | native) is NOT part of the frozen 7-key return -- per
the WU-7 ticket's own "Format" field, leafblower is the only adapter with a
build axis, so `run_leafblower(problem, method, build)` takes it as an
explicit third argument; it surfaces in output only via the weights_ref
filename (runs_schema.json: `build` is a driver-added column, not part of
an adapter's own return schema).

Single-thread BLAS is forced before numpy import (CLAUDE.md determinism
rule -- required for reproducible timing/parity across this repo).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import json  # noqa: E402
import resource  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Optional  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import leafblower  # noqa: E402  -- the ONLY leafblower import: public package entry point

_THIS_DIR = Path(__file__).resolve().parent           # benchmarks/study/python
_STUDY_DIR = _THIS_DIR.parent                           # benchmarks/study
_REPO_ROOT = _STUDY_DIR.parent.parent                    # repo root
sys.path.insert(0, str(_STUDY_DIR / "common"))
from metrics import margin_stats  # noqa: E402

WEIGHTS_DIR = _REPO_ROOT / "weights"

with open(_STUDY_DIR / "spec" / "status_enum.json") as _f:
    _STATUS_ENUM = set(json.load(_f)["$defs"]["StatusEnum"]["enum"])

LBW_METHODS = [
    "oris", "oris_soft", "raking", "sinkhorn", "greenkhorn",
    "chebyshev", "greg", "logit", "newton_kl",
]

# RK_* status code (leafblower.h) -> harmonized enum. Direct 1:1 map for the
# codes harvest() RETURNS normally (0,1,4,5). Codes 2 (infeasible) and 3
# (bad_arg) are ONLY ever raised (RuntimeError / ValueError respectively,
# python/leafblower/_harvest.py:658-668) -- they never reach this map.
_LBW_STATUS_MAP = {0: "converged", 1: "no_conv", 4: "budget", 5: "stall"}


def _peak_rss_bytes() -> int:
    """Best-effort in-process high-water-mark RSS (bytes).

    Linux ru_maxrss is kilobytes (contract.md Section 2.7 describes the
    *subprocess*-level VmHWM sampled by the driver at process exit; this
    in-process ru_maxrss is the adapter-level best-effort proxy, matching
    the convention already established by benchmarks/study/python/competitors.py).
    """
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024)


def _weights_path(solver_id: str, problem_id: str, thread: int, build: str) -> Path:
    WEIGHTS_DIR.mkdir(parents=True, exist_ok=True)
    return WEIGHTS_DIR / f"{solver_id}__{problem_id}__t{thread}__{build}.parquet"


def _write_weights(weights: np.ndarray, solver_id: str, problem_id: str, thread: int, build: str) -> str:
    path = _weights_path(solver_id, problem_id, thread, build)
    pd.DataFrame({"weight": np.asarray(weights, dtype=np.float64)}).to_parquet(path)
    return str(path.relative_to(_REPO_ROOT))


def _nan_sentinel(n: int, solver_id: str, problem_id: str, thread: int, build: str) -> str:
    """contract.md Section 2.1: hard-failure runs still write a length-n
    all-NaN sentinel vector so weights_ref is never dangling."""
    return _write_weights(np.full(n, np.nan), solver_id, problem_id, thread, build)


def _bound_violation(weights: np.ndarray, problem: dict, atol: float = 1e-9) -> bool:
    lo, hi = problem["bounds"]["min"], problem["bounds"]["max"]
    return bool(np.any(weights < lo - atol) or np.any(weights > hi + atol))


def _recompute_converged(weights: np.ndarray, problem: dict) -> bool:
    """contract.md Section 2.4: converged is ALWAYS harness-recomputed here
    via common/metrics.margin_stats' uniform margin-L-infinity check against
    problem['tol'] -- never the solver's own self-report."""
    groups = {m: problem["data"][m].to_numpy() for m in problem["margins"]}
    ms = margin_stats(weights, groups, problem["targets"])
    return bool(ms["margin_linf"] <= problem["tol"])


def _thread() -> int:
    return int(os.environ.get("OMP_NUM_THREADS", "1"))


def run_leafblower(problem: dict, method: str, build: str, thread: Optional[int] = None) -> dict[str, Any]:
    """Run one leafblower method on `problem` through the public
    leafblower.harvest() API and return the frozen 7-key contract v2 result.

    Parameters
    ----------
    problem : problem dict as returned by common/problem_io.py's loader
    method  : one of LBW_METHODS
    build   : "portable" | "native" -- tags the weights_ref filename only
    thread  : thread-count tag for the weights_ref filename; defaults to
              OMP_NUM_THREADS (single-thread-BLAS convention, CLAUDE.md)
    """
    assert method in LBW_METHODS, f"unknown leafblower method {method!r}"
    if thread is None:
        thread = _thread()
    solver_id = f"leafblower_{method}"
    n = len(problem["data"])
    t0 = time.perf_counter()

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        try:
            res = leafblower.harvest(
                data=problem["data"],
                targets=problem["targets"],
                method=method,
                min_weight=problem["bounds"]["min"],
                max_weight=problem["bounds"]["max"],
                design_weights=problem["design_weights"],
                convergence={"absolute": problem["tol"]},
                attach_weights=False,
            )
        except RuntimeError as e:
            # status==2 (infeasible) -- python/leafblower/_harvest.py:658-664
            wall = time.perf_counter() - t0
            ref = _nan_sentinel(n, solver_id, problem["id"], thread, build)
            return dict(weights_ref=ref, iterations=None, status="infeasible",
                        converged=False, error_message=str(e),
                        wall_time_s=float(wall), peak_rss_bytes=_peak_rss_bytes())
        except ValueError as e:
            # status==3 (bad_arg) AND all pre-validation errors -- _harvest.py:667-668
            wall = time.perf_counter() - t0
            ref = _nan_sentinel(n, solver_id, problem["id"], thread, build)
            return dict(weights_ref=ref, iterations=None, status="bad_arg",
                        converged=False, error_message=str(e),
                        wall_time_s=float(wall), peak_rss_bytes=_peak_rss_bytes())
        except Exception as e:  # noqa: BLE001 -- run() must never raise (contract.md Section 1)
            wall = time.perf_counter() - t0
            ref = _nan_sentinel(n, solver_id, problem["id"], thread, build)
            return dict(weights_ref=ref, iterations=None, status="error",
                        converged=False, error_message=f"{type(e).__name__}: {e}",
                        wall_time_s=float(wall), peak_rss_bytes=_peak_rss_bytes())

        result = res["result"]
        status_code = result["status"]
        status = _LBW_STATUS_MAP.get(status_code, "error")  # unmapped code -- defensive only

        weights = np.asarray(res["weights"], dtype=np.float64)
        if _bound_violation(weights, problem):
            status = "bound_violation"

        converged = _recompute_converged(weights, problem)
        wall = time.perf_counter() - t0
        ref = _write_weights(weights, solver_id, problem["id"], thread, build)
        iterations = result.get("iterations")

        msgs = [str(w.message) for w in caught]
        return dict(
            weights_ref=ref,
            iterations=(None if iterations is None else int(iterations)),
            status=status, converged=converged,
            error_message=("; ".join(msgs) if msgs else None),
            wall_time_s=float(wall), peak_rss_bytes=_peak_rss_bytes(),
        )


LEAFBLOWER_PY_ADAPTERS = {
    f"leafblower_{m}_py": (lambda problem, build, _m=m: run_leafblower(problem, _m, build))
    for m in LBW_METHODS
}
