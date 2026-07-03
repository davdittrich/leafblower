"""CR-E12 (leafblower-5ye4.12): the float64 dtype guard in _bindings.cpp was dead
code — py::array::forcecast coerces any numeric weights input to float64 before the
guard ran. Removing it must not change behavior: non-float64 numeric weights are
still accepted (coerced), and the ndim/size guard still rejects wrong shapes."""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import pytest

from leafblower._leafblower import calibrate

RK_ALG_RAKING = 3
RK_OK = 0


def _inputs(n=100, weights=None):
    a = pd.Categorical(np.random.RandomState(1).choice(["x", "y"], n),
                       categories=["x", "y"]).codes.astype(np.int32)
    w = np.ones(n) if weights is None else weights
    return n, [np.ascontiguousarray(a)], [2], [np.array([0.5, 0.5])], w


def _params():
    return {"algorithm": RK_ALG_RAKING, "min_weight": 0.1, "max_weight": 5.0,
            "inner_max_iter": 200, "outer_max_iter": 200, "tol_abs": 1e-6}


def test_int_weights_coerced_via_forcecast():
    # int32 weights (not float64) are forcecast to float64 — accepted, not rejected
    # by the removed dead dtype guard.
    n, gid, cc, tg, _ = _inputs()
    w_int = np.ones(n, dtype=np.int32)
    _s, w_out, r = calibrate(n, 1, w_int, gid, cc, tg, _params(), None)
    assert r["status"] == RK_OK
    assert w_out.sum() > 0


def test_wrong_length_weights_still_rejected():
    n, gid, cc, tg, _ = _inputs()
    with pytest.raises(Exception):
        calibrate(n, 1, np.ones(n - 1), gid, cc, tg, _params(), None)
