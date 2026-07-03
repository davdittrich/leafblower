"""CR-C15b (leafblower-kxna.21): the low-level path (direct c_api / raw calibrate)
bypasses the R/Python wrapper max_iterations>=1 guards. logit with inner_max_iter=0
made kMaxNewtonIters=min(50,0)=0, the Newton loop never ran, and w_best stayed
value-initialized -> degenerate all-zero weights + STALL. The core now rejects it
with BADARG (consistent with the wrapper semantics and chebyshev CR-C16)."""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd

from leafblower._leafblower import calibrate

RK_ALG_LOGIT = 10
RK_ERR_BADARG = 3
RK_OK = 0


def _inputs(n=200):
    a = pd.Categorical(np.random.RandomState(1).choice(["x", "y"], n),
                       categories=["x", "y"]).codes.astype(np.int32)
    return n, [np.ascontiguousarray(a)], [2], [np.array([0.5, 0.5])], np.ones(n)


def test_logit_rejects_zero_inner_budget_not_zero_weights():
    n, gid, cc, tg, w = _inputs()
    p = {"algorithm": RK_ALG_LOGIT, "min_weight": 0.1, "max_weight": 5.0,
         "inner_max_iter": 0, "outer_max_iter": 0, "tol_abs": 1e-6}
    _s, w_out, r = calibrate(n, 1, w.copy(), gid, cc, tg, p, None)
    assert r["status"] == RK_ERR_BADARG
    assert "inner_max_iter" in r["message"]
    # weights are left untouched (input copy), NOT zeroed
    assert w_out.sum() > 0


def test_logit_positive_inner_budget_still_converges():
    n, gid, cc, tg, w = _inputs()
    p = {"algorithm": RK_ALG_LOGIT, "min_weight": 0.1, "max_weight": 5.0,
         "inner_max_iter": 50, "outer_max_iter": 50, "tol_abs": 1e-6}
    _s, _w, r = calibrate(n, 1, w.copy(), gid, cc, tg, p, None)
    assert r["status"] == RK_OK
