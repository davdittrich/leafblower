"""CR-F7b (leafblower-5ye4.16): design_effect must enforce the documented weight
contract (no NA, all finite, sum > 0) on BOTH the 1-arg and 4-arg paths, mirroring
R design_effect.R (dtkn.7). Previously neither Python path validated weights."""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
import pytest

from leafblower import design_effect


def _valid4():
    n = 20
    w = np.ones(n)
    outcome = np.arange(n, dtype=np.float64)
    data = pd.DataFrame({"g": (["x", "y"] * (n // 2))})
    target = {"g": {"x": 0.5, "y": 0.5}}
    return w, outcome, data, target


# ── 1-arg path (Kish deff_K) ─────────────────────────────────────────────────
def test_1arg_nan_rejected():
    w = np.ones(10)
    w[3] = np.nan
    with pytest.raises(ValueError, match="NA"):
        design_effect(w)


def test_1arg_inf_rejected():
    w = np.ones(10)
    w[3] = np.inf
    with pytest.raises(ValueError, match="finite"):
        design_effect(w)


def test_1arg_zero_sum_rejected():
    with pytest.raises(ValueError, match="positive"):
        design_effect(np.zeros(10))


def test_1arg_negative_sum_rejected():
    with pytest.raises(ValueError, match="positive"):
        design_effect(np.full(10, -1.0))


def test_1arg_valid_ok():
    assert abs(design_effect(np.ones(10)) - 1.0) < 1e-12  # equal weights → deff == 1


# ── 5ye4.18: non-numeric weights rejected (parity with R is.numeric) ──
def test_1arg_string_weights_rejected():
    # numpy would silently parse "1.5" → 1.5; R's is.numeric() rejects a character
    # vector, so Python must too.
    with pytest.raises(ValueError, match="numeric"):
        design_effect(["1.0", "2.0", "3.0"])


def test_1arg_bool_weights_rejected():
    # is.numeric(logical) == FALSE in R; np.number excludes bool.
    with pytest.raises(ValueError, match="numeric"):
        design_effect(np.array([True, False, True]))


# ── 4-arg path (Henry & Valliant deff_H) — same contract, was fully bypassed ─
def test_4arg_nan_rejected():
    w, outcome, data, target = _valid4()
    w[2] = np.nan
    with pytest.raises(ValueError, match="NA"):
        design_effect(w, outcome, data, target)


def test_4arg_inf_rejected():
    w, outcome, data, target = _valid4()
    w[2] = np.inf
    with pytest.raises(ValueError, match="finite"):
        design_effect(w, outcome, data, target)


def test_4arg_zero_sum_rejected():
    w, outcome, data, target = _valid4()
    w[:] = 0.0
    with pytest.raises(ValueError, match="positive"):
        design_effect(w, outcome, data, target)


def test_4arg_valid_ok():
    w, outcome, data, target = _valid4()
    d = design_effect(w, outcome, data, target)
    assert np.isfinite(d) and d > 0
