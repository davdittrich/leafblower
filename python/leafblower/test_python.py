import pytest
import numpy as np

def test_harvest_returns_copy():
    """weights_out must be a copy, not a view into input."""
    from leafblower._leafblower import calibrate
    n = 100
    weights = np.ones(n, dtype=np.float64)
    gids = [np.zeros(n, dtype=np.int32)]
    cats = [2]
    tgts = [np.array([0.5, 0.5])]
    # Half in cat 0, half in cat 1
    gids[0][50:] = 1
    status, weights_out, res = calibrate(n, 1, weights, gids, cats, tgts)
    weights_out[0] = 9999.0
    assert weights[0] != 9999.0, "weights_out must be a copy"

def test_convergence_unknown_key_raises():
    from leafblower import harvest
    import pandas as pd
    df = pd.DataFrame({"x": ["a","b","a","b"]})
    tgts = {"x": {"a": 0.5, "b": 0.5}}
    with pytest.raises(ValueError, match="unknown convergence key"):
        harvest(df, tgts, convergence={"bogus_key": 0.01})

def test_min_weight_badarg_python():
    from leafblower import harvest
    import pandas as pd
    df = pd.DataFrame({"x": ["a","b","a","b"]})
    tgts = {"x": {"a": 0.5, "b": 0.5}}
    with pytest.raises(Exception):
        harvest(df, tgts, min_weight=5.0, max_weight=5.0)


def _make_fixture(n=1000):
    """Balanced 2-level fixture for parity tests."""
    import pandas as pd
    half = n // 2
    df = pd.DataFrame({"x": ["a"] * half + ["b"] * half})
    tgts = {"x": {"a": 0.5, "b": 0.5}}
    return df, tgts


def test_default_convergence_is_max_err_improvement():
    """Default convergence: metric=max_err, rule=improvement (WU-G)."""
    from leafblower import harvest
    df, tgts = _make_fixture(n=500)
    res = harvest(df, tgts, max_weight=5, method="ieppa")
    r = res.attrs.get("result", {})
    cu = r.get("convergence_used", {})
    assert cu.get("metric") == "max_err", f"expected max_err, got {cu.get('metric')}"
    assert cu.get("rule") == "improvement", f"expected improvement, got {cu.get('rule')}"


def test_all_metrics_in_result_pct_change_removed():
    """l1_weight_change and grake_norm present; pct_change absent (WU-G)."""
    from leafblower import harvest
    df, tgts = _make_fixture(n=500)
    res = harvest(df, tgts, max_weight=5, method="ieppa")
    r = res.attrs.get("result", {})
    assert "l1_weight_change" in r, "result must contain l1_weight_change"
    assert "grake_norm" in r, "result must contain grake_norm"
    assert "pct_change" not in r, "pct_change must be absent (renamed to l1_weight_change)"


def test_all_6_metrics_present():
    """All 6 quality metrics must be present in result attrs (WU-G)."""
    from leafblower import harvest
    df, tgts = _make_fixture(n=1000)
    res = harvest(df, tgts, max_weight=5, method="ieppa")
    r = res.attrs.get("result", {})
    for key in ("max_error", "mean_error", "kl", "chi2", "l1_weight_change", "grake_norm"):
        assert key in r, f"Missing metric: {key}"


def test_best_error_le_max_error():
    """best_error <= max_error always (best iterate tracked throughout)."""
    from leafblower import harvest
    df, tgts = _make_fixture(n=1000)
    res = harvest(df, tgts, max_weight=5, method="ieppa", max_iterations=200)
    r = res.attrs.get("result", {})
    assert "best_error" in r, "result must contain best_error"
    assert r["best_error"] <= r["max_error"] + 1e-12, (
        f"best_error={r['best_error']} > max_error={r['max_error']}"
    )
