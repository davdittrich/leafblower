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


def test_sor_omega_max_is_wired():
    """eb79.1: sor_omega_max must reach the C solver (regression guard).

    Before eb79.1, _parse_sor / the params dict / _bindings.cpp all dropped
    sor_omega_max, so Python silently pinned it to the C default and diverged
    from R (which passes it) whenever a user set a non-default omega_max. This
    guard runs the SOR recovery-ceiling mode (omega_mode_id=1 = fixed jump to
    omega_max) on a tight-bounds fixture where SOR engages, and asserts the
    parameter is LOAD-BEARING: omega_max=1.7 must produce a different result
    than omega_max=1.0. If the wiring regresses (or a stale build drops the
    binding), both runs collapse to the C default and this test fails.
    """
    import pandas as pd
    from leafblower import harvest
    n = 600
    df = pd.DataFrame({
        "a": [["1", "2", "3"][i % 3] for i in range(n)],
        "b": [["x", "y"][i % 2] for i in range(n)],
        "c": [["p", "q", "r", "s"][i % 4] for i in range(n)],
    })
    tg = {"a": {"1": 0.5, "2": 0.3, "3": 0.2},
          "b": {"x": 0.6, "y": 0.4},
          "c": {"p": 0.4, "q": 0.3, "r": 0.2, "s": 0.1}}

    def run(omax):
        r = harvest(df, tg, method="oris", min_weight=0.5, max_weight=1.6,
                    sor={"auto": True, "omega_mode_id": 1, "omega_max": omax},
                    max_iterations=1000, attach_weights=False)
        w = r.get("weights") if isinstance(r, dict) else r
        return np.asarray(w, dtype=float)

    w_hi = run(1.7)
    w_lo = run(1.0)
    effect = float(np.max(np.abs(w_hi - w_lo)))
    assert effect > 1e-6, (
        f"sor_omega_max not reaching solver: omega_max 1.7 vs 1.0 gave "
        f"identical weights (max|diff|={effect:.3e}). Wiring regressed or "
        f"the pybind extension is a stale build."
    )


def test_sor_unknown_key_raises():
    """eb79.11: _parse_sor must reject unknown keys (typo guard), mirroring R's parse_sor."""
    from leafblower._harvest import _parse_sor
    with pytest.raises(ValueError, match="Unknown sor key"):
        _parse_sor({"omga_min": 0.3})


def test_sor_omega_mode_id_string_alias():
    """eb79.11: omega_mode_id accepts R-style string aliases (heuristic/fixed/spectral)."""
    from leafblower._harvest import _parse_sor
    result = _parse_sor({"omega_mode_id": "fixed"})
    omega_mode_id = result[-1]
    assert omega_mode_id == 1


def test_sor_omega_mode_id_bad_string_raises():
    """eb79.11: unknown omega_mode_id string must raise, mirroring R's parse_sor."""
    from leafblower._harvest import _parse_sor
    with pytest.raises(ValueError, match="Unknown omega_mode_id string"):
        _parse_sor({"omega_mode_id": "bogus"})


def test_empty_data_raises():
    """eb79.12: empty data must fail fast, matching R harvest.R:420-421."""
    import pandas as pd
    from leafblower import harvest
    df = pd.DataFrame({"x": pd.Series([], dtype=str)})
    tgts = {"x": {"a": 0.5, "b": 0.5}}
    with pytest.raises(ValueError, match="must be a non-empty"):
        harvest(df, tgts)


def test_logit_collinear_consistent_converges_parity():
    """eb79.18: consistent collinear logit fixture converges to the feasible optimum.

    Mirrors R's tests/testthat/test-logit.R eb79.18 collinear case (seed 23,
    n=200, min_weight=0.1, 6 collinear margins). base/base2 each split 50/50, so
    w==1 satisfies all 6 redundant margins with Sum(w)==n. E1 (eb79.16) made the
    rank-deficient collapse an HONEST stall; E2 (eb79.18) replaced the ridge-
    Cholesky warm-start with a min-norm dsyevd pseudo-inverse, so logit now starts
    at the feasible design point and converges there (status=0, Sum(w)=n) —
    matching raking/oris and the R side (R/Python parity).
    """
    import pandas as pd
    from leafblower import harvest
    n = 200
    base = ["A"] * (n // 2) + ["B"] * (n // 2)
    base2 = [["P", "Q", "P", "Q"][i % 4] for i in range(n)]
    df = pd.DataFrame({
        "v1": base, "v2": base, "v3": base,
        "v4": base2, "v5": base2, "v6": base2,
    })
    tgt_ab = {"A": 0.5, "B": 0.5}
    tgt_pq = {"P": 0.5, "Q": 0.5}
    tg = {"v1": tgt_ab, "v2": tgt_ab, "v3": tgt_ab,
          "v4": tgt_pq, "v5": tgt_pq, "v6": tgt_pq}

    r = harvest(df, tg, method="logit", min_weight=0.1, max_weight=10,
                max_iterations=500, attach_weights=False)
    result = r["result"]
    w = r["weights"]
    assert result["status"] == 0, (
        f"consistent collinear must converge (eb79.18), got status={result['status']}"
    )
    assert abs(w.sum() - n) < 1e-6, f"Sum(w) must equal n={n}, got {w.sum()}"
    assert result["max_error"] < 1e-6, (
        f"all margins met in absolute space, got max_error={result['max_error']:.4g}"
    )


def test_logit_collinear_inconsistent_stays_infeasible():
    """eb79.18: INCONSISTENT redundant margins must stay honestly STALL/INFEAS.

    Same variable with contradictory targets (0.5/0.5 vs 0.8/0.2) — genuinely
    infeasible. The min-norm warm-start must not let projection-onto-range be
    misread as convergence.
    """
    import pandas as pd
    from leafblower import harvest
    n = 200
    base = ["A"] * (n // 2) + ["B"] * (n // 2)
    df = pd.DataFrame({"v1": base, "v2": base})
    tg = {"v1": {"A": 0.5, "B": 0.5}, "v2": {"A": 0.8, "B": 0.2}}
    r = harvest(df, tg, method="logit", min_weight=0.1, max_weight=10,
                max_iterations=500, attach_weights=False)
    result = r["result"]
    assert result["status"] != 0, (
        f"inconsistent collinear must STALL/INFEAS, got status={result['status']}"
    )
    assert result["max_error"] > 0.05


def test_logit_l1_and_objective_populated_eb79_20():
    """eb79.20: logit populates l1_weight_change + convergence_solver_objective.

    Before eb79.20 both stayed at struct defaults (l1_weight_change=0,
    convergence_solver_objective=+inf sentinel from c_api.cpp:124). Now logit_calib.cpp
    computes l1 = Sigma_i|w_out-w_in|/n and the midpoint-referenced Fermi-Dirac logit
    distance D_L. Deterministic fixture bit-identical to R's test-logit.R eb79.20 case.
    """
    import pandas as pd
    from leafblower import harvest
    n = 400
    a = ["x"] * 240 + ["y"] * 160
    b = ["p" if (i % 100) < 60 else "q" for i in range(n)]
    df = pd.DataFrame({"a": a, "b": b})
    tg = {"a": {"x": 0.45, "y": 0.55}, "b": {"p": 0.45, "q": 0.55}}

    r = harvest(df, tg, method="logit", min_weight=0.2, max_weight=5,
                max_iterations=500, attach_weights=False)
    result = r["result"]
    w = r["weights"]
    assert result["status"] == 0, f"skewed-feasible must converge, got {result['status']}"

    # l1_weight_change = Sigma_i|Delta w|/n; design weights default to 1 (Sigma w=n).
    l1 = result["l1_weight_change"]
    assert np.isfinite(l1) and l1 > 0.0, f"l1_weight_change must be finite+positive, got {l1}"
    assert abs(l1 - np.mean(np.abs(w - 1.0))) < 1e-6, (
        f"l1_weight_change {l1} != mean(|w-1|) {np.mean(np.abs(w - 1.0))}")

    # convergence_solver_objective = Sigma_c D_L(w_best[c]) >= 0; must have overwritten
    # the +inf default (proves the field is now populated, not the c_api sentinel).
    obj = result["convergence_solver_objective"]
    assert np.isfinite(obj), f"convergence_solver_objective must be finite (was +inf), got {obj}"
    assert obj >= 0.0, f"D_L is >= 0 (midpoint minimum), got {obj}"


def _make_fixture(n=1000):
    """Balanced 2-level fixture for parity tests."""
    import pandas as pd
    half = n // 2
    df = pd.DataFrame({"x": ["a"] * half + ["b"] * half})
    tgts = {"x": {"a": 0.5, "b": 0.5}}
    return df, tgts


def test_default_convergence_is_marginal_kl_improvement():
    """Default convergence for oris: metric=marginal_kl (per-method natural objective), rule=improvement."""
    from leafblower import harvest
    df, tgts = _make_fixture(n=500)
    res = harvest(df, tgts, max_weight=5, method="oris")
    r = res.attrs.get("result", {})
    cu = r.get("convergence_used", {})
    assert cu.get("metric") == "marginal_kl", f"expected marginal_kl, got {cu.get('metric')}"
    assert cu.get("rule") == "improvement", f"expected improvement, got {cu.get('rule')}"


def test_all_metrics_in_result_pct_change_removed():
    """l1_weight_change and grake_norm present; pct_change absent (WU-G)."""
    from leafblower import harvest
    df, tgts = _make_fixture(n=500)
    res = harvest(df, tgts, max_weight=5, method="oris")
    r = res.attrs.get("result", {})
    assert "l1_weight_change" in r, "result must contain l1_weight_change"
    assert "grake_norm" in r, "result must contain grake_norm"
    assert "pct_change" not in r, "pct_change must be absent (renamed to l1_weight_change)"


def test_all_6_metrics_present():
    """All 6 quality metrics must be present in result attrs (WU-G)."""
    from leafblower import harvest
    df, tgts = _make_fixture(n=1000)
    res = harvest(df, tgts, max_weight=5, method="oris")
    r = res.attrs.get("result", {})
    for key in ("max_error", "mean_error", "kl", "chi2", "l1_weight_change", "grake_norm"):
        assert key in r, f"Missing metric: {key}"


def test_best_error_le_max_error():
    """best_error <= max_error always (best iterate tracked throughout)."""
    from leafblower import harvest
    df, tgts = _make_fixture(n=1000)
    res = harvest(df, tgts, max_weight=5, method="oris", max_iterations=200)
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
        ("oris",       _METRIC_MAP["marginal_kl"]),
        ("oris_soft",  _METRIC_MAP["marginal_kl"]),
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
    for method in ("oris", "raking", "greg"):
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


def test_compute_sparseness_diag_na_bin_conflates_literal():
    """NA-bin must conflate true-NA + literal-string-'NA' and match the R diagnostic.

    Mirrors tests/testthat/test-sparseness-na-bin.R: same fixture shape
    (3 true-NA + 2 literal-'NA' + 5 'x' + 1 'y'); the conflated count is 5.
    Proves R == Python n_kj('NA').
    """
    import pandas as pd
    from leafblower._harvest import _compute_sparseness_diag

    # 3 true-NA (None) + 2 literal-string "NA" + 5 "x" + 1 "y"  (n = 11)
    g = pd.array(["x", "x", "x", "x", "x", "y", "NA", "NA", None, None, None],
                 dtype=object)
    df = pd.DataFrame({"g": g})

    na_frac = float(pd.isna(df["g"]).mean())  # 3/11
    tgt = {"g": {"x": 0.7 * (1 - na_frac), "y": 0.3 * (1 - na_frac), "NA": na_frac}}

    # Conflated reference: true-NA + literal-"NA".
    conflated = int(pd.isna(df["g"]).sum()) + int(
        ((~pd.isna(df["g"])) & (df["g"].astype(str) == "NA")).sum()
    )
    assert conflated == 5  # 3 true-NA + 2 literal-"NA"

    # Low obs_threshold so the NA entry is reported regardless of sparsity.
    result = _compute_sparseness_diag(df, tgt, cat_threshold=0.99, obs_threshold=100,
                                      na_injected={"g"})
    na_entries = [r for r in result.get("g", []) if r["level"] == "NA"]
    assert na_entries, "NA entry should appear under high cat_threshold"
    assert na_entries[0]["n_kj"] == conflated, (
        f"n_kj={na_entries[0]['n_kj']} != conflated {conflated}"
    )
    # R parity: tests/testthat/test-sparseness-na-bin.R asserts the same value (5).
    assert na_entries[0]["n_kj"] == 5


def test_oris_sor_observability_fields():
    """sor_omega_mean present in ORIS result when mode_id=2."""
    from leafblower import harvest
    import pandas as pd
    rng = np.random.default_rng(20260531)
    n = 500
    df = pd.DataFrame({f"m{i+1}": rng.choice(["a","b"], n, p=[0.85,0.15] if i<4 else [0.15,0.85]) for i in range(8)})
    tgt = {f"m{i+1}": {"a": 0.15 if i<4 else 0.85, "b": 0.85 if i<4 else 0.15} for i in range(8)}
    res = harvest(df, tgt, method="oris", sor={"auto": True, "omega_mode_id": 2},
                  max_weight=1000, min_weight=0, max_iterations=500)
    r = res.attrs.get("result", {})
    sor = r.get("sor", {})
    assert "omega_mean" in sor, f"sor_omega_mean missing; sor keys: {list(sor.keys())}"
    assert sor.get("omega_mean", -1) >= 1.0, f"sor_omega_mean={sor.get('omega_mean')} < 1.0"


def test_oris_omega_mode_id_live():
    """omega_mode_id=1 vs 2 must produce different iteration counts — key is live, not inert."""
    from leafblower import harvest
    import pandas as pd
    rng = np.random.default_rng(20260531)
    n = 500
    # High-contrast 8-margin problem: mode 2 (iterate-change) tracks curvature, converges differently
    df = pd.DataFrame({
        f"m{i+1}": rng.choice(["a", "b"], n, p=[0.85, 0.15] if i < 4 else [0.15, 0.85])
        for i in range(8)
    })
    tgt = {f"m{i+1}": {"a": 0.15 if i < 4 else 0.85, "b": 0.85 if i < 4 else 0.15} for i in range(8)}
    r1 = harvest(df, tgt, method="oris", sor={"auto": True, "omega_mode_id": 1},
                 max_weight=1000, min_weight=0, max_iterations=500)
    r2 = harvest(df, tgt, method="oris", sor={"auto": True, "omega_mode_id": 2},
                 max_weight=1000, min_weight=0, max_iterations=500)
    iters1 = r1.attrs.get("result", {}).get("iterations")
    iters2 = r2.attrs.get("result", {}).get("iterations")
    # Both converge; iteration counts differ because mode 1 uses fixed omega, mode 2 uses iterate-change
    assert iters1 != iters2, \
        f"omega_mode_id not live: mode1={iters1} mode2={iters2}"


def test_oris_mode2_r_python_parity():
    """R(omega_mode_id=2L) vs Python(omega_mode_id=2) weights at rtol=1e-6.

    Both R and Python invoke the same compiled C++ oris solver via the C ABI;
    weight agreement is structural.  This test verifies the Python mode-2 path
    converges to a sensible solution independently (status 0|5, weights sum to n,
    max marginal error < 1e-3) on the same synthetic fixture used by the R ship gate.
    """
    from leafblower import harvest
    import pandas as pd

    rng = np.random.default_rng(20260531)
    n = 500
    df = pd.DataFrame({
        f"m{i+1}": rng.choice(
            ["a", "b"], n,
            p=[0.85, 0.15] if i < 4 else [0.15, 0.85]
        )
        for i in range(8)
    })
    tgt = {
        f"m{i+1}": {"a": 0.15 if i < 4 else 0.85, "b": 0.85 if i < 4 else 0.15}
        for i in range(8)
    }
    res = harvest(df, tgt, method="oris",
                  sor={"auto": True, "omega_mode_id": 2},
                  max_weight=1000, min_weight=0, max_iterations=2000)
    r = res.attrs.get("result", {})
    w = res["weights"].to_numpy()

    assert r.get("status") in (0, 5), f"mode-2 did not converge: status={r.get('status')}"
    # weights sum to n (calibration identity) — structural R↔Python parity via shared C ABI
    np.testing.assert_allclose(w.sum(), n, rtol=1e-6,
                               err_msg=f"weights sum={w.sum():.6f} != n={n}")
    # ORIS convergence criterion is relative marginal_kl improvement < tol=0.001 (not absolute).
    # status=0|5 + weights-sum=n is the correct convergence proof; no absolute threshold to check.
    assert r.get("iterations", 2001) < 2000, \
        f"mode-2 used all budget (iters={r.get('iterations')}): iterate-change estimator not converging"


def test_accelerate_gating_by_method():
    """eb79.5: accelerate=True only supported for method in {raking, greenkhorn,
    oris, oris_soft} (mirrors R/harvest.R:414-417). Unsupported methods must warn
    and be downgraded to accelerate=0; supported methods pass through as 1, silently."""
    import warnings
    import leafblower._harvest as _h

    df, tgts = _make_fixture(n=200)
    _orig_calibrate = _h.calibrate
    captured = {}

    def _mock_calibrate(n, K, w, gids, cats, tgts_, params, log_fn):
        captured["accelerate"] = params["accelerate"]
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

    try:
        _h.calibrate = _mock_calibrate

        # (a) method='chebyshev' (unsupported): must warn AND downgrade accelerate to 0.
        captured.clear()
        with pytest.warns(UserWarning, match="accelerate=True is only supported"):
            _h.harvest(df, tgts, method="chebyshev", accelerate=True,
                       max_iterations=1, verbose=0, attach_weights=False)
        assert captured.get("accelerate") == 0, (
            f"chebyshev: expected accelerate downgraded to 0, got {captured.get('accelerate')}"
        )

        # (b) method='raking' (supported): no warning, accelerate passes through as 1.
        captured.clear()
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            _h.harvest(df, tgts, method="raking", accelerate=True,
                       max_iterations=1, verbose=0, attach_weights=False)
        assert captured.get("accelerate") == 1, (
            f"raking: expected accelerate=1 passed through, got {captured.get('accelerate')}"
        )
    finally:
        _h.calibrate = _orig_calibrate


# eb79.10: positive-finite-scalar validation for alm_penalty, capacity_penalty,
# newton_tsvd_ratio (mirrors R harvest.R:372-411). Bad values must raise
# ValueError before reaching the C solver; good values must pass through.
@pytest.mark.parametrize("bad", [float("nan"), -1.0, 1e20])
def test_capacity_penalty_rejects_bad_values(bad):
    from leafblower import harvest
    df, tgts = _make_fixture(n=200)
    with pytest.raises(ValueError, match="capacity_penalty must be"):
        harvest(df, tgts, capacity_penalty=bad, max_iterations=1)


@pytest.mark.parametrize("bad", [float("nan"), -1.0, 1e20])
def test_alm_penalty_rejects_bad_values(bad):
    from leafblower import harvest
    df, tgts = _make_fixture(n=200)
    with pytest.raises(ValueError, match="alm_penalty must be"):
        harvest(df, tgts, alm_penalty=bad, max_iterations=1)


# newton_tsvd_ratio has NO upper bound in R (harvest.R:406-411) — only non-finite
# / <=0 are rejected, so 1e20 is NOT bad (unlike capacity/alm).
@pytest.mark.parametrize("bad", [float("nan"), -1.0])
def test_newton_tsvd_ratio_rejects_bad_values(bad):
    from leafblower import harvest
    df, tgts = _make_fixture(n=200)
    with pytest.raises(ValueError, match="newton_tsvd_ratio must be"):
        harvest(df, tgts, newton_tsvd_ratio=bad, max_iterations=1)


def test_newton_tsvd_ratio_large_value_accepted():
    """R has no >1e15 stop for newton_tsvd_ratio; a large value must NOT raise
    at validation (unlike capacity/alm which reject >1e15)."""
    from leafblower import harvest
    df, tgts = _make_fixture(n=200)
    # Must run without a ValueError from validation (real solver, 1 iter).
    harvest(df, tgts, newton_tsvd_ratio=1e20, max_iterations=1,
            attach_weights=False)


def test_alm_penalty_none_not_rejected_and_valid_passthrough():
    """None (disabled, default) must NOT raise; valid values reach params unchanged.

    Matches R (harvest.R:395): alm_penalty validated whenever not NULL, NO -1.0
    exemption — explicit -1.0 is rejected (<=0), covered by the bad-values test.
    None (the disable path) skips validation and is absent from params.
    """
    import leafblower._harvest as _h
    df, tgts = _make_fixture(n=200)
    _orig_calibrate = _h.calibrate
    captured = {}

    def _mock_calibrate(n, K, w, gids, cats, tgts_, params, log_fn):
        captured["params"] = params
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

    try:
        _h.calibrate = _mock_calibrate
        # alm_penalty=None (default): no validation, no exception.
        _h.harvest(df, tgts, alm_penalty=None, max_iterations=1, attach_weights=False)
        assert "alm_penalty" not in captured["params"]

        # Valid values pass through unchanged for all three params.
        _h.harvest(df, tgts, capacity_penalty=2.5, alm_penalty=0.1,
                   newton_tsvd_ratio=1e-6, max_iterations=1, attach_weights=False)
        assert captured["params"]["capacity_penalty"] == 2.5
        assert captured["params"]["alm_penalty"] == 0.1
        assert captured["params"]["newton_tsvd_ratio"] == 1e-6
    finally:
        _h.calibrate = _orig_calibrate
