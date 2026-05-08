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
    status, weights_out, res = calibrate(n, 1, weights, gids, cats, tgts, {}, None)
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


def test_default_convergence_is_marginal_kl_improvement():
    """Default convergence for ieppa: metric=marginal_kl (per-method natural objective), rule=improvement."""
    from leafblower import harvest
    df, tgts = _make_fixture(n=500)
    res = harvest(df, tgts, max_weight=5, method="ieppa")
    r = res.attrs.get("result", {})
    cu = r.get("convergence_used", {})
    assert cu.get("metric") == "marginal_kl", f"expected marginal_kl, got {cu.get('metric')}"
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


def test_per_method_default_metric():
    """Per-method natural convergence metric selected when user gives no explicit metric."""
    import pandas as pd
    import numpy as np
    from leafblower._harvest import _METRIC_MAP

    n = 200
    rng = np.random.default_rng(42)
    df = pd.DataFrame({
        "age": rng.choice(["18-34", "35-54", "55+"], n),
        "sex": rng.choice(["M", "F"], n),
    })
    tgt = {
        "age": {"18-34": 0.35, "35-54": 0.40, "55+": 0.25},
        "sex": {"M": 0.48, "F": 0.52},
    }

    cases = [
        ("ieppa",       _METRIC_MAP["marginal_kl"]),
        ("ieppa_soft",  _METRIC_MAP["marginal_kl"]),
        ("raking",      _METRIC_MAP["kl"]),
        ("greenkhorn",  _METRIC_MAP["kl"]),
        ("sinkhorn",    _METRIC_MAP["kl"]),
        ("newton_kl",   _METRIC_MAP["kl"]),
        ("greg",        _METRIC_MAP["chi2"]),
        ("chebyshev",   _METRIC_MAP["max_err"]),
        ("logit",       _METRIC_MAP["max_err"]),
    ]

    import leafblower._harvest as _h
    _orig_calibrate = _h.calibrate
    captured = {}

    def _mock_calibrate(n, K, w, gids, cats, tgts, params, log_fn):
        captured["metric"] = params["metric"]
        weights = np.ones(n, dtype=np.float64)
        result = {
            "status": 0, "iterations": 1, "max_error": 0.0,
            "algorithm_used": 0, "message": "", "n_bounds_violated": 0,
            "n_bounds_clamped": 0, "mean_error": 0.0, "kl": 0.0, "chi2": 0.0,
            "l1_weight_change": 0.0, "grake_norm": 0.0,
            "convergence_metric": params["metric"], "convergence_rule": params["rule"],
            "convergence_tol": params["pct_tol"], "convergence_iter": 1,
            "convergence_solver_objective": 0.0,
            "convergence_minimized_metric": params["metric"],
            "best_error": 0.0, "best_iter": 1,
            "sor_min_omega": 1.0, "sor_n_damped": 0,
            "alm_capacity_mu_final": 0.0, "alm_n_growth_events": 0,
            "alm_max_dual_norm": 0.0, "alm_sum_drift": 0.0,
            "metric_first_check": 0.0, "metric_prev_check": 0.0, "prev_check_iter": 0,
        }
        return 0, weights, result

    for method, expected_metric in cases:
        captured.clear()
        try:
            _h.calibrate = _mock_calibrate  # inside try so restore fires on KI/SIG too
            _h.harvest(df, tgt, method=method, convergence={"tol": 1e-4},
                       max_iterations=1, verbose=0, attach_weights=False)
        except (TypeError, KeyError, IndexError) as exc:
            # post-calibrate result unpacking may fail with mock result — acceptable
            if not captured:
                raise AssertionError(f"mock never called for method={method!r}") from exc
        finally:
            _h.calibrate = _orig_calibrate
        assert captured.get("metric") == expected_metric, (
            f"method={method!r}: expected metric={expected_metric}, got {captured.get('metric')}"
        )

    # Guard: user-explicit metric must NOT be overridden by per-method default
    for method in ("ieppa", "raking", "greg"):
        captured.clear()
        try:
            _h.calibrate = _mock_calibrate
            _h.harvest(df, tgt, method=method,
                       convergence={"metric": "kl", "tol": 1e-4},
                       max_iterations=1, verbose=0, attach_weights=False)
        except (TypeError, KeyError, IndexError):
            pass
        finally:
            _h.calibrate = _orig_calibrate
        assert captured.get("metric") == _METRIC_MAP["kl"], (
            f"explicit metric='kl' overridden for method={method!r}: got {captured.get('metric')}"
        )


def test_compute_sparseness_diag_na_bin():
    """NA-bin level must not be falsely reported as sparse when add_na_proportion=True."""
    import pandas as pd
    import numpy as np
    from leafblower._harvest import _compute_sparseness_diag

    n = 200
    rng = np.random.default_rng(99)
    # 40% NaN — well above obs_threshold=30
    mask = rng.random(n) < 0.4
    strs = rng.choice(["A", "B"], n)
    v1 = pd.array([None if m else s for m, s in zip(mask, strs)], dtype=object)
    df = pd.DataFrame({"v1": v1, "v2": rng.choice(["X", "Y"], n)})

    # Simulate what add_na_proportion=True injects
    na_frac = float(pd.isna(df["v1"]).mean())
    tgt = {
        "v1": {"A": (1 - na_frac) * 0.5, "B": (1 - na_frac) * 0.5, "NA": na_frac},
        "v2": {"X": 0.5, "Y": 0.5},
    }

    na_count = int(pd.isna(df["v1"]).sum())  # ~80 NaN rows
    result = _compute_sparseness_diag(df, tgt, cat_threshold=0.01, obs_threshold=30,
                                      na_injected={"v1"})

    # NA-bin has ~80 obs — must NOT be sparse
    v1_sparse_levels = [r["level"] for r in result.get("v1", [])]
    assert "NA" not in v1_sparse_levels, (
        f"NA-bin falsely flagged sparse: {result.get('v1', [])}"
    )

    # Validate n_kj: use low threshold to force the NA entry to appear, check count
    result_strict = _compute_sparseness_diag(df, tgt, cat_threshold=0.99, obs_threshold=1,
                                             na_injected={"v1"})
    na_entries = [r for r in result_strict.get("v1", []) if r["level"] == "NA"]
    assert na_entries, "NA entry should appear under low cat_threshold"
    assert na_entries[0]["n_kj"] == na_count, (
        f"n_kj={na_entries[0]['n_kj']} != actual NaN count {na_count}"
    )

    # Without na_injected, literal "NA" key uses value_counts (returns 0 for NaN column)
    result_ungated = _compute_sparseness_diag(df, tgt, cat_threshold=0.01, obs_threshold=30)
    ungated_na = [r for r in result_ungated.get("v1", []) if r["level"] == "NA"]
    assert ungated_na, "Without na_injected, NA-bin should be flagged sparse (n_kj=0)"
    assert ungated_na[0]["n_kj"] == 0, "Without na_injected, NA n_kj must be 0 (counts miss NaN)"
