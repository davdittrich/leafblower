"""CR-F9b (leafblower-5ye4.15): _parse_convergence must reject a plateau-rule tol
outside (0,1) on BOTH the explicit `tol` and the `pct` shorthand entry paths,
matching R harvest.R (dtkn.9). Previously Python silently accepted them."""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import warnings

import numpy as np
import pandas as pd
import pytest

from leafblower import harvest
from leafblower._harvest import _parse_convergence


def _fixture(n=100):
    df = pd.DataFrame({"a": np.random.RandomState(1).choice(["x", "y"], n)})
    return df, {"a": {"x": 0.5, "y": 0.5}}, n


def _run(conv):
    df, tg, _ = _fixture()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return harvest(df, tg, convergence=conv, attach_weights=False)


# ── pct shorthand path (the dtkn.9 gap) ──────────────────────────────────────
def test_pct_shorthand_out_of_range_rejected():
    with pytest.raises(ValueError, match=r"pct.*\(0,1\).*plateau"):
        _parse_convergence({"pct": 5})


def test_pct_shorthand_zero_rejected():
    with pytest.raises(ValueError, match=r"\(0,1\)"):
        _parse_convergence({"pct": 0})


def test_pct_shorthand_one_rejected():
    with pytest.raises(ValueError, match=r"\(0,1\)"):
        _parse_convergence({"pct": 1})


def test_pct_shorthand_in_range_accepted():
    pct_tol, _abs, _m, rule, _sw = _parse_convergence({"pct": 0.001})
    assert rule == 2  # plateau
    assert abs(pct_tol - 0.001) < 1e-15


# ── explicit tol path (already checked in R long form) ───────────────────────
def test_explicit_tol_plateau_out_of_range_rejected():
    with pytest.raises(ValueError, match=r"tol.*\(0,1\).*plateau"):
        _parse_convergence({"tol": 5, "rule": "plateau"})


def test_explicit_tol_plateau_in_range_accepted():
    pct_tol, _abs, _m, rule, _sw = _parse_convergence({"tol": 0.01, "rule": "plateau"})
    assert rule == 2
    assert abs(pct_tol - 0.01) < 1e-15


# ── non-plateau rules are unaffected (tol range is a plateau-only contract) ───
def test_threshold_rule_large_tol_unaffected():
    # rule='threshold' with tol=5 is a valid absolute stopping criterion, not plateau.
    _pct, abs_tol, _m, rule, _sw = _parse_convergence({"tol": 5, "rule": "threshold"})
    assert rule == 0
    assert abs(abs_tol - 5.0) < 1e-15


def test_pct_with_explicit_nonplateau_rule_unaffected():
    # pct shorthand but rule overridden to improvement → plateau gate off, no check.
    _pct, _abs, _m, rule, _sw = _parse_convergence({"pct": 5, "rule": "improvement"})
    assert rule == 1


# ── end-to-end through harvest() (both paths surface the error to the user) ───
def test_harvest_pct_out_of_range_raises():
    with pytest.raises(ValueError, match=r"\(0,1\)"):
        _run({"pct": 5})


def test_harvest_pct_in_range_ok():
    out = _run({"pct": 0.001})
    assert out["weights"].shape[0] == 100
