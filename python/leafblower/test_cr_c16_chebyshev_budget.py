"""CR-C16 (leafblower-kxna.16): chebyshev rejects a non-positive inner iteration
budget up front (BADARG + populated message) instead of running the IPM loop zero
times and returning the interior-shifted init with best_error=inf / NOCONV / empty
message. Reached via the low-level c_api binding (the public harvest wrapper already
guards max_iterations>=1)."""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd

from leafblower._leafblower import calibrate

RK_ALG_CHEBYSHEV = 5
RK_ERR_BADARG = 3
RK_OK = 0


def _inputs(n=200):
    a = pd.Categorical(np.random.RandomState(1).choice(["x", "y"], n),
                       categories=["x", "y"]).codes.astype(np.int32)
    return n, [np.ascontiguousarray(a)], [2], [np.array([0.5, 0.5])], np.ones(n)


def test_chebyshev_rejects_zero_inner_budget():
    n, gid, cc, tg, w = _inputs()
    p = {"algorithm": RK_ALG_CHEBYSHEV, "min_weight": 0.1, "max_weight": 5.0,
         "inner_max_iter": 0, "outer_max_iter": 0, "tol_abs": 1e-6}
    status, _w, r = calibrate(n, 1, w.copy(), gid, cc, tg, p, None)
    assert r["status"] == RK_ERR_BADARG
    assert "inner_max_iter" in r["message"]


def test_chebyshev_positive_inner_budget_still_converges():
    n, gid, cc, tg, w = _inputs()
    p = {"algorithm": RK_ALG_CHEBYSHEV, "min_weight": 0.1, "max_weight": 5.0,
         "inner_max_iter": 50, "outer_max_iter": 50, "tol_abs": 1e-6}
    status, _w, r = calibrate(n, 1, w.copy(), gid, cc, tg, p, None)
    assert r["status"] == RK_OK
