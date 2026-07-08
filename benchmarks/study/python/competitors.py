"""benchmarks/study/python/competitors.py -- WU-6 (leafblower-2ouc.7).

Uniform-contract Python adapters wrapping COMPETITOR calibration / raking /
IPF / OT / balancing packages: ipfn, weightipy, POT (ot.sinkhorn/greenkhorn,
K=2 redundant-IPF cross-check), ott-jax (K=2, excluded from ranked timing
per DESIGN.md Section 2), svy (verified maintained successor to archived
samplics -- Gap E), scipy (trust-constr), cvxpy (CLARABEL L-inf baseline),
balance (ipfn-backed rake).

STRICT SEPARATION (user constraint, 2026-07-08): this file imports NOTHING
from leafblower's own src/, r_bridge.cpp, R/, python/leafblower/. It wraps
competitor packages only; leafblower's own adapters are WU-7.

Contract v2 (benchmarks/study/spec/contract.md): every `run(problem)` here
returns EXACTLY
    {weights_ref, iterations, status, converged, error_message,
     wall_time_s, peak_rss_bytes}
`status` is the *package's own* outcome, mapped to the 8-value harmonized
enum (status_enum.json). `converged` is ALWAYS independently recomputed
here from common/metrics.py's margin_stats() against problem["tol"] --
never the package's self-report (contract.md Section 2.4).

Hyperparameters are read verbatim from spec/hyperparams.json's
"Python_competitors" table (frozen, WU-9T) -- never problem["tol"].
problem["tol"] is used ONLY for the harness-side `converged` recomputation.
spec/tol_mapping.json is leafblower-solver-only (its own top-level key is
"leafblower_solvers"); it does not apply to competitor adapters.

Single-thread BLAS is forced before numpy import (CLAUDE.md determinism
rule -- required for reproducible timing/parity across this repo).
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import contextlib  # noqa: E402
import functools  # noqa: E402
import io  # noqa: E402
import json  # noqa: E402
import math  # noqa: E402
import resource  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any, Callable  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

_THIS_DIR = Path(__file__).resolve().parent          # benchmarks/study/python
_STUDY_DIR = _THIS_DIR.parent                          # benchmarks/study
_REPO_ROOT = _STUDY_DIR.parent.parent                   # repo root
sys.path.insert(0, str(_STUDY_DIR / "common"))
from metrics import margin_stats  # noqa: E402

WEIGHTS_DIR = _REPO_ROOT / "weights"

with open(_STUDY_DIR / "spec" / "hyperparams.json") as _f:
    # Each entry in "Python_competitors" is {package, version, arm,
    # entrypoint, params, convergence_bounds_quantity, citation, note} --
    # adapters below read solver *knobs* only, so _HP flattens straight to
    # each solver's "params" sub-dict.
    _HP = {k: v["params"] for k, v in json.load(_f)["Python_competitors"].items()}

with open(_STUDY_DIR / "spec" / "status_enum.json") as _f:
    _STATUS_ENUM = set(json.load(_f)["$defs"]["StatusEnum"]["enum"])

# Pinned OT entropic epsilon (WU-6 decision -- no separate "OT-reduction
# spec" doc exists yet; recorded here per DESIGN.md Section 2's "adapter
# golden must pin this" directive). This value is a MATHEMATICAL identity,
# not a tuning knob: entropic-regularized OT with cost c_ij = -log(T_ij)
# (T = the observed 2-margin seed cell-count table) and regularization eps
# has kernel K_ij = exp(-c_ij/eps) = T_ij**(1/eps). Sinkhorn's fixed point
# on kernel K is exactly IPF's fixed point on K's cell table (Sinkhorn IS
# IPF, applied multiplicatively via potentials u,v instead of direct
# row/col rescaling). That equivalence K == T -- and hence Sinkhorn's
# coupling P_ij = u_i K_ij v_j reproducing IPF's converged table exactly --
# holds ONLY at eps == 1 (K_ij = T_ij**1). eps != 1 raises T to a power and
# silently diverges from IPF (empirically confirmed: eps=0.5 and eps=2.0
# both produce >0.4 max-abs-error per-cell-multiplier vs IPF on a 2x2
# interacting golden, vs ~3e-16 at eps=1.0). Structural zero cells (T_ij=0,
# c_ij=-log(0)=+inf in principle) are clamped to _OT_BIG_COST for numerical
# finiteness. Validated against a hand-solved 2-margin IPF golden in
# test_competitors.py::test_ot_extraction_matches_ipf_golden.
_OT_EPSILON = 1.0
_OT_BIG_COST = 700.0  # exp(-700) underflows to 0.0 in float64; finite and safe


# ---------------------------------------------------------------------------
# Shared plumbing
# ---------------------------------------------------------------------------

def _peak_rss_bytes() -> int:
    """Best-effort in-process high-water-mark RSS.

    Linux ru_maxrss is kilobytes (contract.md Section 2.7 describes the
    *subprocess*-level VmHWM sampled by the driver at process exit; this
    in-process ru_maxrss is the adapter-level best-effort proxy used when
    `run()` is not itself the sole content of its subprocess).
    """
    return int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss * 1024)


def _weights_path(solver_id: str, problem_id: str, thread: int = 1, build: str = "na") -> Path:
    WEIGHTS_DIR.mkdir(parents=True, exist_ok=True)
    return WEIGHTS_DIR / f"{solver_id}__{problem_id}__t{thread}__{build}.parquet"


def _write_weights(weights: np.ndarray, solver_id: str, problem_id: str) -> str:
    path = _weights_path(solver_id, problem_id)
    pd.DataFrame({"weight": np.asarray(weights, dtype=np.float64)}).to_parquet(path)
    return str(path.relative_to(_REPO_ROOT))


def _nan_sentinel(n: int, solver_id: str, problem_id: str) -> str:
    """contract.md Section 2.1: hard-failure runs still write a length-n
    all-NaN sentinel vector so weights_ref is never dangling."""
    return _write_weights(np.full(n, np.nan), solver_id, problem_id)


def _recompute_converged(weights: np.ndarray, problem: dict) -> tuple[bool, float]:
    """contract.md Section 2.4: converged is ALWAYS harness-recomputed here
    via common/metrics.margin_stats' uniform margin-L-infinity check against
    problem['tol'] -- never the package's own self-report."""
    groups = {m: problem["data"][m].to_numpy() for m in problem["margins"]}
    ms = margin_stats(weights, groups, problem["targets"])
    return bool(ms["margin_linf"] <= problem["tol"]), float(ms["margin_linf"])


def _bound_violation(weights: np.ndarray, problem: dict, atol: float = 1e-9) -> bool:
    lo, hi = problem["bounds"]["min"], problem["bounds"]["max"]
    return bool(np.any(weights < lo - atol) or np.any(weights > hi + atol))


def _margin_constraint_matrix(problem: dict) -> tuple[sp.csr_matrix, np.ndarray]:
    """Sparse linear equality system A w = rhs encoding every margin's
    every category (count-scale: rhs_j = T_kj * W, W = sum(design_weights)).
    Used by the two NLP/QP-based adapters (scipy trust-constr, cvxpy).
    """
    df = problem["data"]
    W = float(problem["design_weights"].sum())
    rows, rhs = [], []
    for m in problem["margins"]:
        cats = df[m].to_numpy()
        for lv, frac in problem["targets"][m].items():
            rows.append(sp.csr_matrix((cats == lv).astype(np.float64)))
            rhs.append(frac * W)
    A = sp.vstack(rows, format="csr")
    return A, np.asarray(rhs, dtype=np.float64)


def _ot_seed_cost(cats1: np.ndarray, cats2: np.ndarray, levels1: list, levels2: list,
                   design_weights: np.ndarray) -> np.ndarray:
    """2-margin OT reduction cost matrix (DESIGN.md Section 2, corrected):
    c_ij = -log(T_ij) where T_ij is the observed seed cell weight (sum of
    design_weights over rows in cell (i,j)). Combined with _OT_EPSILON==1,
    this makes the Sinkhorn kernel K_ij = exp(-c_ij) == T_ij exactly, so
    Sinkhorn's fixed point reproduces IPF's fixed point on the same seed
    table (see _OT_EPSILON docstring above). Structural zeros (T_ij == 0)
    get _OT_BIG_COST (finite stand-in for -log(0) = +inf)."""
    idx1 = {lv: i for i, lv in enumerate(levels1)}
    idx2 = {lv: i for i, lv in enumerate(levels2)}
    T = np.zeros((len(levels1), len(levels2)), dtype=np.float64)
    for c1, c2, d in zip(cats1, cats2, design_weights):
        if c1 in idx1 and c2 in idx2:
            T[idx1[c1], idx2[c2]] += d
    return np.where(T > 0.0, -np.log(np.maximum(T, 1e-300)), _OT_BIG_COST)


def _ot_weights_from_coupling(u: np.ndarray, v: np.ndarray, cats1: np.ndarray, cats2: np.ndarray,
                               levels1: list, levels2: list, design_weights: np.ndarray) -> np.ndarray:
    """w_i = d_i * u_{cat1(i)} * v_{cat2(i)} -- the FULL Sinkhorn-coupling
    cell multiplier (both potentials), matching IPF's biproportional
    per-cell factor r_i*c_j exactly (empirically confirmed to ~1e-16 vs a
    hand-solved interacting 2x2 golden). A left-potential-only extraction
    (u alone) is WRONG whenever the column margin also needs adjustment --
    it silently drops the v_j factor and only happens to coincide with IPF
    when the column margin is already balanced. Renormalized so
    sum(w) == sum(d) (leafblower's Sigma-w=n convention, CLAUDE.md) to
    remove the (u,v) -> (u*t, v/t) scale ambiguity intrinsic to Sinkhorn
    potentials."""
    idx1 = {lv: i for i, lv in enumerate(levels1)}
    idx2 = {lv: i for i, lv in enumerate(levels2)}
    row = np.array([idx1[c] for c in cats1])
    col = np.array([idx2[c] for c in cats2])
    w = u[row] * v[col] * design_weights
    W = float(design_weights.sum())
    s = float(w.sum())
    if s > 0:
        w = w * (W / s)
    return w


def _adapt(solver_id: str):
    """Decorator: uniform timing/exception/weights-write/converged-recompute
    envelope shared by every adapter below. `solve_fn(problem)` returns
    `(weights_or_None, iterations_or_None, status, error_message_or_None)`;
    weights=None is reserved for hard bad_arg/infeasible failures with no
    usable iterate (writes the NaN sentinel instead)."""

    def decorator(solve_fn: Callable[[dict], tuple]):
        @functools.wraps(solve_fn)
        def run(problem: dict) -> dict[str, Any]:
            t0 = time.perf_counter()
            n = len(problem["data"])
            try:
                weights, iterations, status, error_message = solve_fn(problem)
                if status not in _STATUS_ENUM:
                    raise AssertionError(f"{solver_id}: unmapped status {status!r}")
                if weights is None:
                    wall = time.perf_counter() - t0
                    ref = _nan_sentinel(n, solver_id, problem["id"])
                    return dict(weights_ref=ref, iterations=None, status=status,
                                converged=False, error_message=error_message,
                                wall_time_s=float(wall), peak_rss_bytes=_peak_rss_bytes())
                weights = np.asarray(weights, dtype=np.float64)
                if _bound_violation(weights, problem):
                    status = "bound_violation"
                converged, _linf = _recompute_converged(weights, problem)
                wall = time.perf_counter() - t0
                ref = _write_weights(weights, solver_id, problem["id"])
                return dict(
                    weights_ref=ref,
                    iterations=(None if iterations is None else int(iterations)),
                    status=status, converged=converged, error_message=error_message,
                    wall_time_s=float(wall), peak_rss_bytes=_peak_rss_bytes(),
                )
            except Exception as e:  # noqa: BLE001 -- run() must never raise (contract.md Section 1)
                wall = time.perf_counter() - t0
                ref = _nan_sentinel(n, solver_id, problem["id"])
                return dict(weights_ref=ref, iterations=None, status="error",
                            converged=False, error_message=f"{type(e).__name__}: {e}",
                            wall_time_s=float(wall), peak_rss_bytes=_peak_rss_bytes())

        run.solver_id = solver_id
        return run

    return decorator


# ---------------------------------------------------------------------------
# ipfn (registry id "ipfn"; KL family, unbounded only)
# Home-turf golden: ipfn's own IPFN(df, aggregates, dimensions) DataFrame
# API on a tiny 2-margin table -- README/PyPI worked example shape.
# ---------------------------------------------------------------------------

@_adapt("ipfn")
def run_ipfn(problem: dict) -> tuple:
    from ipfn.ipfn import ipfn as IPFN

    if problem["bounds"]["min"] > 0.0 or math.isfinite(problem["bounds"]["max"]):
        return None, None, "bad_arg", "ipfn is KL/IPF, unbounded-only per registry.json"

    df = problem["data"]
    margins = problem["margins"]
    hp = _HP["ipfn"]

    work = df[margins].copy()
    work["_d"] = problem["design_weights"]
    cell = work.groupby(margins, observed=True)["_d"].sum().reset_index(name="total")
    seed = cell["total"].to_numpy().copy()

    W = float(problem["design_weights"].sum())
    aggregates, dimensions = [], []
    for m in margins:
        aggregates.append(pd.Series({k: v * W for k, v in problem["targets"][m].items()},
                                     name="total", dtype=float))
        dimensions.append([m])

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ipf = IPFN(cell, aggregates, dimensions, weight_col="total",
                   convergence_rate=hp["convergence_rate"],
                   max_iteration=hp["max_iteration"],
                   rate_tolerance=hp["rate_tolerance"], verbose=2)
        cell_out, self_converged, conv_df = ipf.iteration()

    mult = np.where(seed > 0, cell_out["total"].to_numpy() / seed, 1.0)
    lookup = cell[margins].copy()
    lookup["_m"] = mult
    merged = work.merge(lookup, on=margins, how="left")
    weights = merged["_m"].fillna(1.0).to_numpy() * problem["design_weights"]

    status = "converged" if self_converged == 1 else "no_conv"
    return weights, len(conv_df), status, None


# ---------------------------------------------------------------------------
# weightipy (registry id "weightipy"; RIM raking, KL family, unbounded only)
# Home-turf golden: weightipy's own scheme_from_dict + weight_dataframe
# README example (single-margin age/gender rim scheme).
# ---------------------------------------------------------------------------

@_adapt("weightipy")
def run_weightipy(problem: dict) -> tuple:
    import weightipy

    if problem["bounds"]["min"] > 0.0 or math.isfinite(problem["bounds"]["max"]):
        return None, None, "bad_arg", "weightipy RIM is unbounded per registry.json"

    df = problem["data"]
    margins = problem["margins"]
    hp = _HP["weightipy"]

    distributions = {m: {k: v * 100.0 for k, v in problem["targets"][m].items()} for m in margins}
    scheme = weightipy.scheme_from_dict(
        distributions,
        rim_params=dict(max_iterations=hp["max_iterations"], convcrit=hp["convcrit"]),
    )
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        out = weightipy.weight_dataframe(df[margins], scheme, weight_column="weights", verbose=False)

    rim_w = out["weights"].to_numpy()
    # weightipy's RIM has no native design-weight input hook: post-multiply
    # by d_i (WU-6 limitation, documented; RIM itself fit on unweighted counts).
    weights = rim_w * problem["design_weights"]

    group = next(iter(scheme.groups.values()))
    iterations = group.get("iterations")
    status = "no_conv" if (iterations is not None and iterations >= hp["max_iterations"]) else "converged"
    return weights, iterations, status, None


# ---------------------------------------------------------------------------
# POT sinkhorn / greenkhorn (registry ids "pot_sinkhorn"/"pot_greenkhorn";
# OT family, K_max=2 -- redundant-IPF cross-check, DESIGN.md Section 2).
# Home-turf golden: POT's own ot.sinkhorn(a, b, M, reg) two-Gaussian tutorial
# shape (a, b on a small support, dense cost matrix).
# ---------------------------------------------------------------------------

def _pot_run(problem: dict, method: str, hp_key: str, solver_id: str) -> tuple:
    import ot

    margins = problem["margins"]
    if len(margins) != 2:
        return None, None, "bad_arg", (
            f"{solver_id}: OT reduction requires K==2 (registry K_max=2); got K={len(margins)}"
        )
    m1, m2 = margins
    df = problem["data"]
    levels1 = list(problem["targets"][m1].keys())
    levels2 = list(problem["targets"][m2].keys())
    a = np.array([problem["targets"][m1][k] for k in levels1], dtype=np.float64)
    b = np.array([problem["targets"][m2][k] for k in levels2], dtype=np.float64)
    M = _ot_seed_cost(df[m1].to_numpy(), df[m2].to_numpy(), levels1, levels2, problem["design_weights"])

    hp = _HP[hp_key]
    fn = ot.sinkhorn if method == "sinkhorn" else ot.bregman.greenkhorn
    G, log = fn(a, b, M, reg=_OT_EPSILON, numItermax=hp["numItermax"], stopThr=hp["stopThr"], log=True)

    u = np.asarray(log["u"], dtype=np.float64)
    v = np.asarray(log["v"], dtype=np.float64)
    weights = _ot_weights_from_coupling(u, v, df[m1].to_numpy(), df[m2].to_numpy(),
                                         levels1, levels2, problem["design_weights"])
    niter = log.get("niter", log.get("n_iter"))  # sinkhorn: "niter"; greenkhorn: "n_iter"
    status = "no_conv" if (niter is not None and niter >= hp["numItermax"]) else "converged"
    return weights, niter, status, None


@_adapt("pot_sinkhorn")
def run_pot_sinkhorn(problem: dict) -> tuple:
    return _pot_run(problem, "sinkhorn", "pot_sinkhorn", "pot_sinkhorn")


@_adapt("pot_greenkhorn")
def run_pot_greenkhorn(problem: dict) -> tuple:
    return _pot_run(problem, "greenkhorn", "pot_greenkhorn", "pot_greenkhorn")


# ---------------------------------------------------------------------------
# ott-jax sinkhorn (registry id "ott_jax_sinkhorn"; OT family, K_max=2).
# Per DESIGN.md Section 2 / contract DoD: wrapped + home-turf golden'd like
# every other competitor, but EXCLUDED FROM RANKED TIMING by the WU-8 run
# matrix (XLA/JAX thread-pinning under {1,4}-thread OMP/BLAS caps is not
# apples-to-apples with the other arms) -- that exclusion is a run-matrix
# (WU-8) concern; this adapter itself is fully functional.
# Home-turf golden: ott-jax's own pointcloud/Sinkhorn quickstart shape,
# adapted here to the same support-cost 2-margin reduction as POT so both
# OT arms are pinned against the identical IPF golden.
# ---------------------------------------------------------------------------

@_adapt("ott_jax_sinkhorn")
def run_ott_jax_sinkhorn(problem: dict) -> tuple:
    import jax.numpy as jnp
    from ott.geometry import geometry
    from ott.problems.linear import linear_problem
    from ott.solvers.linear import sinkhorn

    margins = problem["margins"]
    if len(margins) != 2:
        return None, None, "bad_arg", (
            f"ott_jax_sinkhorn: OT reduction requires K==2 (registry K_max=2); got K={len(margins)}"
        )
    m1, m2 = margins
    df = problem["data"]
    levels1 = list(problem["targets"][m1].keys())
    levels2 = list(problem["targets"][m2].keys())
    a = jnp.asarray([problem["targets"][m1][k] for k in levels1], dtype=jnp.float64)
    b = jnp.asarray([problem["targets"][m2][k] for k in levels2], dtype=jnp.float64)
    M = _ot_seed_cost(df[m1].to_numpy(), df[m2].to_numpy(), levels1, levels2, problem["design_weights"])

    hp = _HP["ott_jax_sinkhorn"]
    geom = geometry.Geometry(cost_matrix=jnp.asarray(M), epsilon=_OT_EPSILON)
    lp = linear_problem.LinearProblem(geom, a=a, b=b)
    solver = sinkhorn.Sinkhorn(threshold=hp["threshold"], min_iterations=hp["min_iterations"],
                                max_iterations=hp["max_iterations"], inner_iterations=hp["inner_iterations"])
    out = solver(lp)

    u = np.asarray(out.scalings[0], dtype=np.float64)
    v = np.asarray(out.scalings[1], dtype=np.float64)
    weights = _ot_weights_from_coupling(u, v, df[m1].to_numpy(), df[m2].to_numpy(),
                                         levels1, levels2, problem["design_weights"])
    status = "converged" if bool(out.converged) else "no_conv"
    return weights, int(out.n_iters), status, None


# ---------------------------------------------------------------------------
# svy (registry id "svy"; verified maintained successor to archived
# samplics -- Gap E: importing samplics itself now emits
# "samplics is archived and no longer maintained. Migrate to 'svy'").
# KL/chi2 families, unbounded only.
# Home-turf golden: svy's own weighting.raking.rake() tutorial shape
# (Design(wgt=...) + Sample + controls dict of category counts).
# ---------------------------------------------------------------------------

@_adapt("svy")
def run_svy(problem: dict) -> tuple:
    import polars as pl
    from svy import Design, Sample
    from svy.weighting.raking import rake

    if problem["bounds"]["min"] > 0.0 or math.isfinite(problem["bounds"]["max"]):
        return None, None, "bad_arg", "svy.weighting.raking.rake is unbounded per registry.json"

    df = problem["data"]
    margins = problem["margins"]
    hp = _HP["svy"]
    W = float(problem["design_weights"].sum())

    pl_df = pl.from_pandas(df[margins].copy())
    pl_df = pl_df.with_columns(pl.Series("_d", problem["design_weights"]))
    sample = Sample(pl_df, Design(wgt="_d"))
    controls = {m: {k: v * W for k, v in problem["targets"][m].items()} for m in margins}

    out = rake(sample, controls=controls, wgt_name="rk_wgt", tol=hp["rake_tol"], max_iter=100)
    weights = out.data["rk_wgt"].to_numpy().astype(np.float64)
    # svy's public rake() API exposes no iteration count (hyperparams.json
    # flags max_iterations/calibrate_tol as undocumented for this call).
    return weights, None, "converged", None


# ---------------------------------------------------------------------------
# scipy trust-constr (registry id "scipy_trust_constr"; newton_kl family,
# bounded). Objective = KL(w || d); constraints = every margin-category
# equality (count scale); bounds = problem bounds.
# Home-turf golden: scipy.optimize.minimize(method='trust-constr') own
# quickstart (quadratic objective + LinearConstraint + Bounds).
# ---------------------------------------------------------------------------

@_adapt("scipy_trust_constr")
def run_scipy_trust_constr(problem: dict) -> tuple:
    from scipy.optimize import Bounds, LinearConstraint, minimize

    d = problem["design_weights"]
    n = len(d)
    hp = _HP["scipy_trust_constr"]
    eps = 1e-12

    def f(w):
        wc = np.maximum(w, eps)
        return float(np.sum(wc * np.log(wc / d)))

    def grad(w):
        wc = np.maximum(w, eps)
        return np.log(wc / d) + 1.0

    def hess(w):
        wc = np.maximum(w, eps)
        return sp.diags(1.0 / wc)

    A, rhs = _margin_constraint_matrix(problem)
    lc = LinearConstraint(A, rhs, rhs)
    lo = problem["bounds"]["min"]
    hi = problem["bounds"]["max"] if math.isfinite(problem["bounds"]["max"]) else np.inf
    bounds = Bounds(np.full(n, lo), np.full(n, hi))
    x0 = np.clip(d, max(lo, eps), hi if math.isfinite(hi) else None)

    res = minimize(f, x0, jac=grad, hess=hess, constraints=[lc], bounds=bounds,
                    method="trust-constr",
                    options=dict(gtol=hp["gtol"], xtol=hp["xtol"],
                                 barrier_tol=hp["barrier_tol"], maxiter=hp["maxiter"]))

    # scipy trust-constr status: 0=maxiter, 1=gtol satisfied, 2=xtol
    # satisfied, 3=callback-requested stop (AdapterStatusMapping, this WU).
    status_map = {0: "budget", 1: "converged", 2: "converged", 3: "no_conv"}
    status = status_map.get(res.status, "no_conv" if res.success else "error")
    return res.x, int(res.nit), status, (None if res.success else str(res.message))


# ---------------------------------------------------------------------------
# cvxpy CLARABEL L-infinity baseline (registry id "cvxpy_linf"; minimax
# family, bounded). Objective = ||w - d||_inf; constraints = every
# margin-category equality (count scale); bounds = problem bounds.
# Home-turf golden: cvxpy's own norm_inf-minimization LP quickstart.
# ---------------------------------------------------------------------------

@_adapt("cvxpy_linf")
def run_cvxpy_linf(problem: dict) -> tuple:
    import cvxpy as cp

    d = problem["design_weights"]
    n = len(d)
    hp = _HP["cvxpy_linf"]

    w = cp.Variable(n)
    A, rhs = _margin_constraint_matrix(problem)
    lo = problem["bounds"]["min"]
    hi = problem["bounds"]["max"]
    constraints = [A @ w == rhs, w >= lo]
    if math.isfinite(hi):
        constraints.append(w <= hi)
    prob = cp.Problem(cp.Minimize(cp.norm_inf(w - d)), constraints)

    # AdapterStatusMapping (this WU): CVXPY solver status string -> enum.
    status_map = {
        cp.OPTIMAL: "converged", cp.OPTIMAL_INACCURATE: "no_conv",
        cp.INFEASIBLE: "infeasible", cp.INFEASIBLE_INACCURATE: "infeasible",
        cp.UNBOUNDED: "error", cp.UNBOUNDED_INACCURATE: "error",
    }
    prob.solve(solver=cp.CLARABEL,
               tol_gap_abs=hp["tol_gap_abs"], tol_gap_rel=hp["tol_gap_rel"],
               tol_feas=hp["tol_feas"], max_iter=hp["max_iter"])

    status = status_map.get(prob.status, "error")
    if w.value is None:
        return None, None, status, f"cvxpy status={prob.status}"
    iterations = None
    if prob.solver_stats is not None:
        iterations = prob.solver_stats.num_iters
    return np.asarray(w.value, dtype=np.float64), iterations, status, None


# ---------------------------------------------------------------------------
# balance (registry id "balance"; KL family via ipfn backend, unbounded).
# Home-turf golden: balance's own Sample.from_frame + set_target +
# adjust(method='rake') README shape.
# ---------------------------------------------------------------------------

@_adapt("balance")
def run_balance(problem: dict) -> tuple:
    from balance.sample_class import Sample

    if problem["bounds"]["min"] > 0.0 or math.isfinite(problem["bounds"]["max"]):
        return None, None, "bad_arg", "balance rake (ipfn backend) is unbounded per registry.json"

    df = problem["data"]
    margins = problem["margins"]
    n = len(df)
    hp = _HP["balance"]

    sample_df = df[margins].copy()
    sample_df.insert(0, "_id", [str(i) for i in range(n)])
    sample_df["_weight"] = problem["design_weights"]

    # Synthetic target frame: independent per-margin resampling honoring
    # each margin's target proportion (rounded to a fixed target_n), which
    # is exactly the "marginal-only" information balance's IPF-backed rake
    # consumes (no joint-target claim is made or needed).
    target_n = max(1000, n)
    rng_target = {}
    for m in margins:
        levels = list(problem["targets"][m].keys())
        fracs = np.array([problem["targets"][m][k] for k in levels])
        counts = np.round(fracs * target_n).astype(int)
        counts[-1] += target_n - counts.sum()
        rng_target[m] = np.repeat(levels, np.maximum(counts, 0))
        rng_target[m] = np.resize(rng_target[m], target_n)
    target_df = pd.DataFrame(rng_target)
    target_df.insert(0, "_id", [str(i) for i in range(target_n)])
    target_df["_weight"] = 1.0

    kwargs = dict(id_column="_id", weight_column="_weight", covar_columns=list(margins))
    s = Sample.from_frame(sample_df, **kwargs)
    t = Sample.from_frame(target_df, **kwargs)
    s2 = s.set_target(t)

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
        adjusted = s2.adjust(method=hp["method"])
    weights = adjusted.df["_weight"].to_numpy().astype(np.float64)
    # balance's rake() returns weights on the TARGET frame's total-weight
    # scale (here target_n, since target_df._weight==1.0 per row), not the
    # sample's own design-weight scale -- rescale to Sigma-w=n (CLAUDE.md
    # convention, matches every other adapter in this file).
    W = float(problem["design_weights"].sum())
    s = float(weights.sum())
    if s > 0:
        weights = weights * (W / s)
    return weights, None, "converged", None


# ---------------------------------------------------------------------------
ADAPTERS: dict[str, Callable[[dict], dict[str, Any]]] = {
    "ipfn": run_ipfn,
    "weightipy": run_weightipy,
    "pot_sinkhorn": run_pot_sinkhorn,
    "pot_greenkhorn": run_pot_greenkhorn,
    "ott_jax_sinkhorn": run_ott_jax_sinkhorn,
    "svy": run_svy,
    "scipy_trust_constr": run_scipy_trust_constr,
    "cvxpy_linf": run_cvxpy_linf,
    "balance": run_balance,
}
