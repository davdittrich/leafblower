"""CR-D6 (leafblower-j7x8.6): newton_kl populates an identical, complete result-field
set whether reached via explicit RK_ALG_NEWTON_KL or via the AUTO fallback.

Before the fix the explicit c_api branch omitted mean_error/kl/chi2/l1_weight_change/
grake_norm/convergence_solver_objective/convergence_minimized_metric (leaving
convergence_solver_objective at its Inf sentinel), and the AUTO-fallback branch
dropped grake_norm and left the abandoned ORIS primary's sor_*/alm_*/homotopy_*
diagnostics stale. Both c_api paths now route through the shared
pack_newton_result_c helper, so the field set is complete and the ORIS-only fields
read their documented non-ORIS defaults.

The explicit branch is reachable through the public harvest() API. AUTO is R-only at
the wrapper layer, so the c_api AUTO-fallback branch is exercised via the low-level
_leafblower.calibrate binding with algorithm=0 (RK_ALG_AUTO).
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import warnings

import numpy as np
import pandas as pd

from leafblower import harvest
from leafblower._leafblower import calibrate

RK_ALG_NEWTON_KL = 11
RK_ALG_AUTO = 0

# ORIS-only diagnostics and their documented non-ORIS defaults (leafblower.h).
# newton_kl owns none of these, so both the explicit and fallback paths must
# report exactly these values.
ORIS_DEFAULTS = {
    "min_alpha_seen": 1.0,
    "final_alpha": 1.0,
    "homotopy_final_factor": 1.0,
    "homotopy_levels_used": 0,
    "greedy_sweeps_taken": 0,
    "eta_final": 0.0,
    "sor_min_omega": 1.0,
    "sor_omega_mean": 1.0,
    "sor_n_damped": 0,
    "sor_any_latched": 0,
    "alm_capacity_mu_final": 0.0,
    "alm_n_growth_events": 0,
    "alm_max_dual_norm": 0.0,
    "alm_sum_drift": 0.0,
}


def _assert_complete_newton_fields(r, label):
    # Base fields the explicit branch previously omitted must now be populated.
    assert np.isfinite(r["convergence_solver_objective"]), (
        f"{label}: convergence_solver_objective must be finite (was Inf/omitted "
        f"pre-D6), got {r['convergence_solver_objective']}"
    )
    for k in ("mean_error", "kl", "chi2", "l1_weight_change", "grake_norm",
              "convergence_minimized_metric"):
        assert k in r, f"{label}: missing base field {k}"
        assert np.isfinite(r[k]), f"{label}: base field {k} not finite: {r[k]}"
    # ORIS-only diagnostics at documented non-ORIS defaults (no stale ORIS state).
    for k, v in ORIS_DEFAULTS.items():
        assert r[k] == v, f"{label}: ORIS-only {k}={r[k]}, expected default {v}"


def test_explicit_newton_kl_populates_complete_fields():
    rng = np.random.RandomState(1)
    n = 300
    df = pd.DataFrame({
        "a": rng.choice(["x", "y"], n),
        "b": rng.choice(["p", "q"], n),
        "c": rng.choice(["m", "o"], n),
    })
    tg = {"a": {"x": .5, "y": .5}, "b": {"p": .5, "q": .5}, "c": {"m": .5, "o": .5}}
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        r = harvest(df, tg, method="newton_kl", max_iterations=50,
                    attach_weights=False)["result"]
    assert r["algorithm_used"] == RK_ALG_NEWTON_KL
    _assert_complete_newton_fields(r, "explicit")


def test_auto_fallback_newton_kl_via_c_api_zeroes_oris_extras():
    # Skewed compressed problem + accelerate + tiny budget: AUTO routes to an
    # ORIS+SRAA primary (which populates ORIS-only diagnostics), NOCONV/BUDGETs,
    # then falls back to newton_kl. Driven through the raw c_api binding because
    # "auto" is disabled in the harvest() wrapper.
    rng = np.random.RandomState(11)
    n = 500
    cols = {
        "a": rng.choice(["x", "y", "z"], n, p=[.85, .10, .05]),
        "b": rng.choice(["p", "q"], n, p=[.75, .25]),
        "c": rng.choice(["m", "o"], n, p=[.6, .4]),
    }
    targets = {
        "a": [0.34, 0.33, 0.33],
        "b": [0.5, 0.5],
        "c": [0.5, 0.5],
    }
    levels = {"a": ["x", "y", "z"], "b": ["p", "q"], "c": ["m", "o"]}
    group_ids = []
    cat_counts = []
    targets_list = []
    for col, lv in levels.items():
        code = pd.Categorical(cols[col], categories=lv).codes.astype(np.int32)
        group_ids.append(np.ascontiguousarray(code))
        cat_counts.append(len(lv))
        targets_list.append(np.asarray(targets[col], dtype=np.float64))
    w = np.ones(n, dtype=np.float64)
    params = {
        "algorithm": RK_ALG_AUTO,
        "min_weight": 0.0,
        "max_weight": 5.0,
        "inner_max_iter": 8,
        "outer_max_iter": 8,
        "tol_abs": 1e-6,
        "accelerate": 1,
    }
    status, _w_out, r = calibrate(
        n, len(levels), w, group_ids, cat_counts, targets_list, params, None
    )
    assert r["algorithm_used"] == RK_ALG_NEWTON_KL, "AUTO must fall back to newton_kl"
    # Same completeness contract as the explicit branch, and no stale ORIS state
    # (homotopy_levels_used would be >=1 if the ORIS primary's value leaked through).
    _assert_complete_newton_fields(r, "auto-fallback")
