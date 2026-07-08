# benchmarks/study/common/ref_convex.py
#
# Independent high-precision convex reference solver -- the RQ5 agreement
# anchor (DESIGN.md Section 6 / Blocker G). Ticket leafblower-2ouc.5 (WU-4).
# Python mirror of ref_convex.R -- see that file's header for the full
# per-family objective derivation and citations against leafblower's own
# C++ solver source (src/raking.cpp, src/logit_calib.cpp, src/chebyshev.cpp).
#
# INDEPENDENCE: this module calls cvxpy only. It NEVER imports leafblower or
# calls harvest() -- it is the correctness anchor, not a benchmarked
# competitor.
#
# Scope limit: a 1e-12 convex solve is infeasible at stepstone scale
# (n ~ 1.58M) -- REF_MAX_N below refuses any larger problem loudly rather
# than faking a result (DESIGN.md Section 6: "1.58M has no anchor").

from __future__ import annotations

import os

# Single-thread BLAS forced BEFORE numpy import for R<->Python parity
# determinism (project convention; see CLAUDE.md), matching metrics.py.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import json  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import cvxpy as cp  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from metrics import margin_stats  # noqa: E402

REF_SOLVER_TOL = 1e-12
REF_MAX_N = 50_000
REF_STRICTLY_CONVEX_FAMILIES = ("kl", "chi2", "logit")
REF_FAMILIES = ("kl", "chi2", "logit", "minimax")
REF_ANCHOR_KIND = {
    "kl": "weight_vector",
    "chi2": "weight_vector",
    "logit": "weight_vector",
    "minimax": "objective_value",
}


def _assert_scope(problem: dict[str, Any]) -> None:
    n = len(problem["data"])
    if n > REF_MAX_N:
        raise ValueError(
            f"ref_convex: problem '{problem['id']}' has n={n} > REF_MAX_N="
            f"{REF_MAX_N} -- a ~1e-12 convex reference solve is infeasible "
            "at this scale (DESIGN.md Section 6: '1.58M has no anchor', "
            "stated limitation, never faked). Exclude from the RQ5 anchor."
        )


def _margin_targets(problem: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Population total N and per-margin category target COUNTS (T_kj*N),
    plus the per-observation category label vector for each margin. Shared
    by every family: the margin constraint is always
    Sum_{i in cat kj} w_i == T_kj * N, linear regardless of family."""
    N = float(np.sum(problem["design_weights"]))
    out = {}
    for m in problem["margins"]:
        groups = problem["data"][m].to_numpy()
        target = problem["targets"][m]
        out[m] = dict(
            groups=groups,
            target_counts={lvl: p * N for lvl, p in target.items()},
        )
    return out


def _objective(family: str, w: cp.Variable, d: np.ndarray, L: float, U: float):
    if family == "kl":
        if np.any(d <= 0):
            raise ValueError(
                "ref_convex: kl family requires d_i > 0 for all i (got "
                "d_i<=0); exclude zero-design-weight rows before calling "
                "solve_ref()."
            )
        return cp.sum(cp.kl_div(w, d))
    elif family == "chi2":
        if np.any(d <= 0):
            raise ValueError(
                "ref_convex: chi2 family requires d_i > 0 for all i (got "
                "d_i<=0); exclude zero-design-weight rows before calling "
                "solve_ref()."
            )
        return cp.sum(cp.multiply(cp.square(w - d), 1.0 / d))
    elif family == "logit":
        if not np.isfinite(U):
            raise ValueError(
                "ref_convex: logit family requires a finite max bound "
                "(matches leafblower's own guard, src/logit_calib.cpp:46 "
                "'max_weight must be finite and positive')."
            )
        return -cp.sum(cp.entr(w - L)) - cp.sum(cp.entr(U - w))
    else:
        raise ValueError(f"ref_convex: unsupported strictly-convex family: {family}")


def solve_ref(problem: dict[str, Any], family: str, tol: float = REF_SOLVER_TOL) -> dict[str, Any]:
    """Solve the independent convex reference problem for one (problem,
    family). Strictly-convex families (kl, chi2, logit) return
    mode="weight_vector" (unique optimum). minimax returns
    mode="objective_value" ONLY (Blocker G: non-unique optimum face) --
    weights are used internally to recompute the achieved L-inf error via
    the shared golden metrics.py::margin_stats(), never stored as an
    anchor."""
    if family not in REF_FAMILIES:
        raise ValueError(
            f"ref_convex: unknown family '{family}' (expected one of: "
            f"{', '.join(REF_FAMILIES)})"
        )
    _assert_scope(problem)

    n = len(problem["data"])
    d = np.asarray(problem["design_weights"], dtype=float)
    L = problem["bounds"]["min"]
    U = problem["bounds"]["max"]
    N = float(np.sum(d))
    mt = _margin_targets(problem)

    w = cp.Variable(n)
    base_cons = [w >= L, cp.sum(w) == N]
    if np.isfinite(U):
        base_cons.append(w <= U)

    eq_cons = []
    for mm in mt.values():
        groups = mm["groups"]
        for lvl, tc in mm["target_counts"].items():
            eq_cons.append(cp.sum(w[groups == lvl]) == tc)

    solve_kwargs = dict(
        solver=cp.CLARABEL,
        tol_gap_abs=tol, tol_gap_rel=tol, tol_feas=tol,
        tol_infeas_abs=tol, tol_infeas_rel=tol,
    )

    if family in REF_STRICTLY_CONVEX_FAMILIES:
        obj = _objective(family, w, d, L, U)
        prob = cp.Problem(cp.Minimize(obj), base_cons + eq_cons)
        prob.solve(**solve_kwargs)
        st = prob.status
        if st not in ("optimal", "optimal_inaccurate"):
            raise RuntimeError(
                f"ref_convex: CLARABEL did not reach optimality (family="
                f"{family}, problem={problem['id']}): status={st}"
            )
        return dict(family=family, mode="weight_vector",
                    weights=np.asarray(w.value, dtype=float),
                    objective=float(prob.value), solver_status=st)
    else:  # minimax
        t = cp.Variable(1)
        linf_cons = []
        for mm in mt.values():
            groups = mm["groups"]
            for lvl, tc in mm["target_counts"].items():
                S = cp.sum(w[groups == lvl])
                linf_cons.append(S - tc <= t * N)
                linf_cons.append(tc - S <= t * N)
        prob = cp.Problem(cp.Minimize(t), base_cons + [t >= 0] + linf_cons)
        prob.solve(**solve_kwargs)
        st = prob.status
        if st not in ("optimal", "optimal_inaccurate"):
            raise RuntimeError(
                "ref_convex: CLARABEL did not reach optimality (family="
                f"minimax, problem={problem['id']}): status={st}"
            )
        w_solved = np.asarray(w.value, dtype=float)
        groups = {m: problem["data"][m].to_numpy() for m in problem["margins"]}
        ms = margin_stats(w_solved, groups, problem["targets"])
        return dict(family="minimax", mode="objective_value",
                    obj_val=ms["margin_linf"], weights=w_solved, solver_status=st)


def ref_weights_path(family: str, problem_id: str, out_dir: str = "benchmarks/study/results") -> Path:
    return Path(out_dir) / "weights" / f"ref_{family}__{problem_id}__t1__na.parquet"


def ref_objective_path(family: str, problem_id: str, out_dir: str = "benchmarks/study/results") -> Path:
    return Path(out_dir) / "ref_objective" / f"{family}__{problem_id}.json"


def store_ref(problem: dict[str, Any], family: str, result: dict[str, Any],
               out_dir: str = "benchmarks/study/results") -> Path:
    """Store a solve_ref() result as a pseudo-solver row (strictly-convex
    families: weights/ parquet, one column `weight`, length n) or the
    achieved-objective-value record (minimax: small JSON -- never a
    weight-vector parquet, so minimax can never be accidentally consumed by
    RQ5 vector correlation)."""
    import pandas as pd

    if result["mode"] == "weight_vector":
        path = ref_weights_path(family, problem["id"], out_dir)
        path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({"weight": result["weights"]}).to_parquet(path, index=False)
        return path
    else:
        path = ref_objective_path(family, problem["id"], out_dir)
        path.parent.mkdir(parents=True, exist_ok=True)
        record = dict(
            problem=problem["id"], family=family, anchor_kind="objective_value",
            obj_val=result["obj_val"], solver_status=result["solver_status"],
            note=("Blocker G: L-inf optimum is a non-unique face; excluded "
                  "from weight-vector correlation (DESIGN.md Section 6)."),
        )
        with open(path, "w") as f:
            json.dump(record, f, indent=2)
        return path
