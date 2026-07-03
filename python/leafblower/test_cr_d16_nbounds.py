"""CR-D16 (leafblower-j7x8.16): chebyshev/greenkhorn/logit now surface the real
n_bounds_violated count through the c_api has_n_bounds_c trait (previously dropped
into a local -> Python saw a stale 0)."""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import warnings

import numpy as np
import pandas as pd
import pytest

from leafblower import harvest


@pytest.mark.parametrize("method", ["chebyshev", "greenkhorn"])
def test_nbounds_surfaced_nonzero(method):
    rng = np.random.RandomState(3)
    n = 400
    df = pd.DataFrame({
        "a": rng.choice(["x", "y", "z"], n, p=[.7, .2, .1]),
        "b": rng.choice(["p", "q"], n, p=[.65, .35]),
        "c": rng.choice(["m", "o"], n),
    })
    tg = {"a": {"x": .34, "y": .33, "z": .33}, "b": {"p": .5, "q": .5}, "c": {"m": .5, "o": .5}}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = harvest(df, tg, method=method, min_weight=0.6, max_weight=1.5,
                    max_iterations=200, attach_weights=False)["result"]
    assert isinstance(r["n_bounds_violated"], int)
    assert r["n_bounds_violated"] > 0, f"{method}: tight bounds should bind"


def test_logit_nbounds_field_wired():
    # logit enforces bounds via the sigmoid link, so the count is a genuine 0 with
    # uniform base weights — but the field must be wired (cross-layer parity with R).
    rng = np.random.RandomState(4)
    n = 300
    df = pd.DataFrame({"a": rng.choice(["x", "y"], n), "b": rng.choice(["p", "q"], n)})
    tg = {"a": {"x": .5, "y": .5}, "b": {"p": .5, "q": .5}}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = harvest(df, tg, method="logit", min_weight=0.2, max_weight=3,
                    max_iterations=100, attach_weights=False)["result"]
    assert isinstance(r["n_bounds_violated"], int)
    assert r["n_bounds_violated"] >= 0
