"""CR-C10b (leafblower-kxna.22): the shared mark_converged helper reported
convergence_tol=pct_tol=0.0 when NEITHER cfg tol was set, even though the run
converged via the st.tol_abs fallback (check_convergence: `curr < tol_abs_fallback`).
mark_converged now reports tol_abs_fallback in that case. Other cases unchanged."""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd

from leafblower._leafblower import calibrate

RK_ALG_RAKING = 3
RK_OK = 0


def _inputs(n=300):
    r1 = pd.Categorical(np.random.RandomState(2).choice(["x", "y"], n),
                        categories=["x", "y"]).codes.astype(np.int32)
    r2 = pd.Categorical(np.random.RandomState(3).choice(["p", "q"], n),
                        categories=["p", "q"]).codes.astype(np.int32)
    gid = [np.ascontiguousarray(r1), np.ascontiguousarray(r2)]
    tg = [np.array([0.5, 0.5]), np.array([0.5, 0.5])]
    return n, gid, [2, 2], tg, np.ones(n)


def test_neither_tol_reports_tol_abs_not_zero():
    n, gid, cc, tg, w = _inputs()
    p = {"algorithm": RK_ALG_RAKING, "min_weight": 0.1, "max_weight": 5.0,
         "inner_max_iter": 500, "outer_max_iter": 500,
         "pct_tol": 0.0, "absolute_tol": 0.0, "tol_abs": 1e-6, "metric": 0, "rule": 0}
    _s, _w, r = calibrate(n, 2, w.copy(), gid, cc, tg, p, None)
    assert r["status"] == RK_OK
    # neither cfg tol set -> report the actual fallback that governed convergence
    assert r["convergence_tol"] == 1e-6, r["convergence_tol"]


def test_absolute_only_still_reports_absolute_tol():
    # Regression guard: the abs-only case must keep reporting absolute_tol (not the
    # fallback), matching test-convergence-criteria.R.
    n, gid, cc, tg, w = _inputs()
    p = {"algorithm": RK_ALG_RAKING, "min_weight": 0.1, "max_weight": 5.0,
         "inner_max_iter": 500, "outer_max_iter": 500,
         "pct_tol": 0.0, "absolute_tol": 1e-5, "tol_abs": 1e-6, "metric": 0, "rule": 0}
    _s, _w, r = calibrate(n, 2, w.copy(), gid, cc, tg, p, None)
    assert r["status"] == RK_OK
    assert r["convergence_tol"] == 1e-5, r["convergence_tol"]
