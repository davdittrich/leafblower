"""CR-E10b (leafblower-5ye4.14): start_weights validation is now identical in R and
Python (canonical contract: reject negatives / non-finite / wrong-length; length-1
broadcasts to n)."""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import warnings

import numpy as np
import pandas as pd
import pytest

from leafblower import harvest


def _fixture(n=100):
    df = pd.DataFrame({"a": np.random.RandomState(1).choice(["x", "y"], n)})
    return df, {"a": {"x": 0.5, "y": 0.5}}, n


def _run(sw):
    df, tg, n = _fixture()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return harvest(df, tg, start_weights=sw, attach_weights=False)


def test_rejects_negative_start_weights():
    _, _, n = _fixture()
    with pytest.raises(ValueError, match="negative"):
        _run(np.array([2.0] * (n - 1) + [-0.5]))


def test_rejects_non_finite_start_weights():
    _, _, n = _fixture()
    with pytest.raises(ValueError, match="non-finite"):
        _run(np.array([1.0] * (n - 1) + [np.nan]))


def test_rejects_wrong_length_start_weights():
    with pytest.raises(ValueError, match="length"):
        _run(np.ones(50))  # neither n nor 1


def test_rejects_tiny_positive_sum():
    _, _, n = _fixture()
    with pytest.raises(ValueError, match="positive value"):
        _run(np.full(n, 1e-18))


def test_length_one_start_weights_broadcasts():
    # length-1 broadcasts to n (parity with R rep()) — must NOT raise
    _run(np.array([3.0]))


def test_bare_scalar_broadcasts():
    # bare Python float broadcasts to n (parity with R scalar) — must NOT raise
    _run(3.0)


def test_rejects_2d_start_weights():
    _, _, n = _fixture()
    with pytest.raises(ValueError, match="1-D"):
        _run(np.ones((n, 1)))


def test_accepts_valid_start_weights():
    _, _, n = _fixture()
    _run(np.ones(n))
