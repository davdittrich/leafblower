"""Surface tests for design_effect Python wrapper."""

import numpy as np
import pandas as pd
import pytest

from leafblower import design_effect, effective_sample_size


def test_design_effect_1arg_matches_kish():
    w = np.array([1.0, 2.0, 3.0, 4.0])
    expected = float(len(w)) * float(np.sum(w**2)) / float(np.sum(w))**2
    assert design_effect(w) == pytest.approx(expected, rel=1e-12)


def test_effective_sample_size():
    w = np.array([1.0, 2.0, 3.0, 4.0])
    assert effective_sample_size(w) == pytest.approx(len(w) / design_effect(w), rel=1e-12)


def test_design_effect_4arg_empty_target_returns_deff_K():
    w = np.array([1.0, 2.0, 3.0, 4.0])
    y = np.array([10.0, 20.0, 30.0, 40.0])
    data = pd.DataFrame({"g": ["a", "b", "a", "b"]})
    d = design_effect(w, outcome=y, data=data, target={})
    assert d == pytest.approx(design_effect(w), rel=1e-12)


def test_design_effect_4arg_r_parity():
    """R-parity: same fixture as test-design.R analytic case -> same deff_H value."""
    np.random.seed(42)
    n = 30
    # 3-level group; level A has code 0 (dropped in C), B=1, C=2.
    groups = np.repeat(["A", "B", "C"], n // 3)
    w = np.ones(n)
    y = np.where(groups == "A", 10.0, np.where(groups == "B", 20.0, 30.0))
    data = pd.DataFrame({"g": groups})
    target = {"g": {"A": 1/3, "B": 1/3, "C": 1/3}}
    d = design_effect(w, outcome=y, data=data, target=target)
    # Analytic: beta_hat = [20, 30] (slopes for B, C vs A); u = [0, 0, 0] (perfect fit).
    # sigma2_u = 0 -> deff_H = 0 is NOT expected here because A group (code 0, dropped)
    # has residual u_A = y_A - 0 = 10 (intercept not included in X).
    # Expected: deff_H approximately 1/3 (analytically derived, same as R test).
    assert d == pytest.approx(1.0/3.0, rel=1e-6)
