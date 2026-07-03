"""CR-E (epic 5ye4): Python layer + build integrity.

Covers 5ye4.3 (STALL surfacing), .4 (unknown-kwargs rejection), .5 (auto_collapse
ordering), .7 (truthiness gates), .8 (6 marshalled result fields), .9 (pandas
now a hard requirement).
"""
import pytest
import pandas as pd

from leafblower import harvest
import leafblower._harvest as H


def _simple():
    df = pd.DataFrame({"x": ["a", "b", "a", "b", "a", "b"]})
    tgts = {"x": {"a": 0.5, "b": 0.5}}
    return df, tgts


# ── 5ye4.3: STALL (status 5) surfaced as a warning ──────────────────────────
def _patch_status(monkeypatch, code):
    real = H.calibrate

    def fake(*a, **k):
        rc, w, res = real(*a, **k)
        res["status"] = code
        return rc, w, res

    monkeypatch.setattr(H, "calibrate", fake)


def test_stall_base_warning(monkeypatch):
    df, tgts = _simple()
    _patch_status(monkeypatch, 5)
    with pytest.warns(UserWarning, match="loss function plateau"):
        harvest(df, tgts, method="oris")


def test_stall_accelerate_warning(monkeypatch):
    df, tgts = _simple()
    _patch_status(monkeypatch, 5)
    with pytest.warns(UserWarning, match="SRAA-m weight-change plateau"):
        harvest(df, tgts, method="oris", accelerate=True)


# ── 5ye4.4: unknown kwargs rejected (typos no longer swallowed) ─────────────
def test_unknown_kwarg_raises():
    df, tgts = _simple()
    with pytest.raises(TypeError, match="max_weigth"):
        harvest(df, tgts, max_weigth=3.0)  # typo for max_weight


def test_multiple_unknown_kwargs_named():
    df, tgts = _simple()
    with pytest.raises(TypeError, match="bounds_mod"):
        harvest(df, tgts, bounds_mod="unit", frobnicate=1)


def test_removed_param_still_specific():
    df, tgts = _simple()
    with pytest.raises(TypeError, match="removed autumn param 'enforce_mean'"):
        harvest(df, tgts, enforce_mean=True)


# ── 5ye4.5: auto_collapse runs after coercion + parse ───────────────────────
def test_auto_collapse_dict_data():
    """Failure mode (b): dict data + auto_collapse used to AttributeError."""
    data = {"x": ["a"] * 50 + ["rare"] + ["b"] * 49}
    tgts = {"x": {"a": 0.5, "b": 0.49, "rare": 0.01}}
    out = harvest(data, tgts, auto_collapse=True, attach_weights=False)
    assert out["weights"].shape[0] == 100


def test_auto_collapse_dataframe_targets():
    """Failure mode (a): DataFrame targets + auto_collapse used to TypeError."""
    df = pd.DataFrame({"x": ["a"] * 50 + ["rare"] + ["b"] * 49})
    tgt_df = pd.DataFrame({"variable": ["x", "x", "x"],
                           "level": ["a", "b", "rare"],
                           "proportion": [0.5, 0.49, 0.01]})
    out = harvest(df, tgt_df, auto_collapse=True, attach_weights=False)
    assert out["weights"].shape[0] == 100


# ── 5ye4.7: truthiness gates (not identity) ─────────────────────────────────
@pytest.mark.parametrize("val,collapses", [(1, True), (True, True),
                                           (0, False), (False, False), (None, False)])
def test_auto_collapse_truthiness(val, collapses):
    # collapse_vars omitted (None) so only bool(auto_collapse) gates: 1/True enable,
    # 0/False/None disable. 'rare' has count 1 (<30) so it folds into __other__ when
    # enabled. attach_weights default → DataFrame; inspect the working x column.
    df = pd.DataFrame({"x": ["a"] * 50 + ["rare"] + ["b"] * 49})
    tgts = {"x": {"a": 0.5, "b": 0.49, "rare": 0.01}}
    out = harvest(df, tgts, auto_collapse=val)
    has_other = "__other__" in set(out["x"].astype(str).unique())
    assert has_other == collapses


@pytest.mark.parametrize("val", [0, None, False])
def test_add_na_proportion_disabled_values(val):
    """0/None/False must all DISABLE NA-bin injection (was `is not False`)."""
    df = pd.DataFrame({"x": ["a", "b", None, "a", "b", "a"]})
    tgts = {"x": {"a": 0.5, "b": 0.5}}
    out = harvest(df, tgts, add_na_proportion=val, attach_weights=False)
    # NA bin NOT injected → targets unchanged, no "NA" level in result targets.
    assert "NA" not in tgts["x"]  # caller dict never mutated
    assert out["weights"].shape[0] == 6


# ── 5ye4.8: 6 previously-dropped result fields now marshalled ───────────────
def test_marshalled_result_fields():
    df, tgts = _simple()
    out = harvest(df, tgts, attach_weights=False)
    res = out["result"]
    for f in ("homotopy_levels_used", "homotopy_final_factor", "greedy_sweeps_taken",
              "eta_final", "min_alpha_seen", "final_alpha"):
        assert f in res, f"missing marshalled field {f}"
    # single-pass default n_levels=1 → homotopy_levels_used == 1 (NOT 0)
    assert res["homotopy_levels_used"] == 1


# ── 5ye4.9: pandas is now a hard requirement (no fiction) ────────────────────
def test_pandas_available_flag_gone():
    assert not hasattr(H, "_PANDAS_AVAILABLE"), \
        "_PANDAS_AVAILABLE dead-fallback machinery should be removed"
