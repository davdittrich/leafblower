#!/usr/bin/env python3
"""benchmarks/study/common/test_metrics.py

Hand-computed golden for benchmarks/study/common/metrics.py (ticket leafblower-2ouc.2).
Same toy and same golden literals as test_metrics.R -- this file is also the
R<->Python parity oracle: every numeric literal here was independently
derived once (see the scratch derivation in the ticket work log) and is
reused verbatim in both languages.

Run: python/.venv/bin/python benchmarks/study/common/test_metrics.py
(single-thread BLAS forced before numpy import, matching metrics.py)
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import math  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
from metrics import (  # noqa: E402
    agreement,
    compute_metrics,
    margin_stats,
)

fail_count = 0


def check(desc, got, want, tol=1e-9):
    if math.isinf(want):
        ok = math.isinf(got) and (got > 0) == (want > 0)
    else:
        ok = math.isclose(got, want, rel_tol=0, abs_tol=tol) or math.isclose(got, want, rel_tol=tol)
    print(f"[{'PASS' if ok else 'FAIL'}] {desc:<60} got={got!r} want={want!r}")
    global fail_count
    if not ok:
        fail_count += 1
    return ok


def check_true(desc, got):
    ok = bool(got)
    print(f"[{'PASS' if ok else 'FAIL'}] {desc:<60} got={got!r}")
    global fail_count
    if not ok:
        fail_count += 1
    return ok


# -----------------------------------------------------------------------------
# Toy: 3x4 -- margin A (3 levels) x margin B (4 levels), n=6 observations.
# Identical to test_metrics.R.
# -----------------------------------------------------------------------------
groupA = np.array(["a1", "a1", "a2", "a2", "a3", "a3"])
groupB = np.array(["b1", "b2", "b3", "b4", "b1", "b2"])
TA = {"a1": 0.3, "a2": 0.3, "a3": 0.4}
TB = {"b1": 0.25, "b2": 0.25, "b3": 0.25, "b4": 0.25}
groups = {"A": groupA, "B": groupB}
targets = {"A": TA, "B": TB}


# ---- Legacy reference impl (verbatim, cited; NOT modified in situ) --------

def legacy_compute_metrics(weights, groupA, groupB, TA, TB):
    """Verbatim algorithm from benchmarks/python_ipf_benchmark.py:31-58
    (compute_metrics), adapted only to take (weights, groupA, groupB, TA, TB)
    instead of (weights, df_pd, tgt) -- same per-category loop and formulas.
    """
    weights = np.asarray(weights, dtype=float)
    W = float(weights.sum())
    n = len(weights)
    max_err = L1 = chi2 = marg_kl = 0.0
    for group, tgt in ((groupA, TA), (groupB, TB)):
        for cat, t in tgt.items():
            S = float(weights[group == cat].sum()) / W
            e = abs(S - t)
            if e > max_err:
                max_err = e
            L1 += e
            exp = t * W
            if exp > 0:
                chi2 += ((S * W - exp) ** 2) / exp
            if t > 0 and S > 0:
                marg_kl += t * np.log(t / S)
    pos = weights > 0
    weight_kl = float(np.sum(weights[pos] * np.log(weights[pos])) / n)  # BUG: hardcodes d_i=1
    return dict(max_err=max_err, L1=L1, chi2=chi2, marg_kl=marg_kl, weight_kl=weight_kl)


# -----------------------------------------------------------------------------
# Case 0: baseline (uniform w, d_i=1) -- sanity + legacy agreement on marg_kl
# -----------------------------------------------------------------------------
print("== Case 0: baseline uniform toy ==")
w0 = np.array([1.0, 1, 1, 1, 1, 1])
m0 = compute_metrics(w0, groups, targets, family="kl")

check("marg_kl_mean (base)", m0["marg_kl_mean"], 0.03430191557553895)
check("marg_kl_max (base)", m0["marg_kl_max"], 0.05889151782819174)
check("margin_linf (base)", m0["margin_linf"], 0.08333333333333334)
check("margin_l1 (base)", m0["margin_l1"], 0.4666666666666667)
check("ESS (base, uniform)", m0["ESS"], 6.0)
check("DEFF_kish (base, uniform)", m0["DEFF_kish"], 1.0)
check("weight_kl (base, d=1, uniform)", m0["weight_kl"], 0.0)
check_true("g_weighted FALSE when d all 1", not m0["g_weighted"])

total_kl_base = 0.009712313322886162 + 0.05889151782819174
check("new impl total (marg_kl_mean*2) matches independent total (base)",
      m0["marg_kl_mean"] * 2, total_kl_base)
lg0 = legacy_compute_metrics(w0, groupA, groupB, TA, TB)
check("legacy python_ipf_benchmark marg_kl matches independent total (base, non-divergent)",
      lg0["marg_kl"], total_kl_base)

# -----------------------------------------------------------------------------
# Case (a): starved category -- assert marg_kl -> +inf, divergent=True, and
# that the legacy `t>0 and S>0` gate SILENTLY DROPS the divergent term (bug).
# -----------------------------------------------------------------------------
print("\n== Case (a): starved category (b4 starved to w=0) ==")
w_starved = np.array([1.0, 1, 1, 0, 1, 2])  # sum(w)=6=n; b4's sole member (obs4) starved
m_starved = compute_metrics(w_starved, groups, targets, family="kl")

check("marg_kl_max is +inf on starved category", m_starved["marg_kl_max"], math.inf)
check("marg_kl_mean is +inf on starved category (propagates, not averaged away)",
      m_starved["marg_kl_mean"], math.inf)
check_true("divergent flag set True", m_starved["marg_kl_divergent"])
check_true("divergent_margins == ['B'] (the starved margin)",
           m_starved["marg_kl_divergent_margins"] == ["B"])

ms_starved = margin_stats(w_starved, groups, targets)
check_true("margin A stays finite (only B starved)", math.isfinite(ms_starved["per_margin"]["A"]["kl"]))
check("margin A KL value (starved case, still finite)",
      ms_starved["per_margin"]["A"]["kl"], 0.05547042424760395)

lg_starved = legacy_compute_metrics(w_starved, groupA, groupB, TA, TB)
check_true("BUG REPRODUCED: legacy python_ipf_benchmark marg_kl is finite (silently drops divergent term)",
           math.isfinite(lg_starved["marg_kl"]))

# -----------------------------------------------------------------------------
# Case (b): d_i != 1 -- assert g-weight DEFF != Kish DEFF, and that legacy
# weight_kl (hardcoded d_i=1) diverges from the new d-aware weight_kl (bug).
# -----------------------------------------------------------------------------
print("\n== Case (b): non-uniform design weights d_i != 1 ==")
w_b = np.array([1.0, 1, 1, 1, 2, 2])
d_b = np.array([1.0, 1, 3, 3, 1, 1])
m_b = compute_metrics(w_b, groups, targets, d=d_b, family="kl")

check("ESS (case b)", m_b["ESS"], 16 / 3)
check("DEFF_kish (case b, raw w)", m_b["DEFF_kish"], 1.125)
check("DEFF_g (case b, g=w/d)", m_b["DEFF_g"], 1.38)
check_true("g-weight DEFF != Kish DEFF (Gap D)", abs(m_b["DEFF_g"] - m_b["DEFF_kish"]) > 1e-6)
check_true("g_weighted True when d_i != 1", m_b["g_weighted"])
check("weight_kl (case b, real d_i)", m_b["weight_kl"], 0.5753641449035616)

lg_b = legacy_compute_metrics(w_b, groupA, groupB, TA, TB)
check("BUG REPRODUCED: legacy python_ipf_benchmark weight_kl (wrongly assumes d_i=1)",
      lg_b["weight_kl"], 0.46209812037329684)
check_true("legacy weight_kl (d_i=1 assumed) != correct d-aware weight_kl",
           abs(lg_b["weight_kl"] - m_b["weight_kl"]) > 1e-3)

bv = compute_metrics(w_b, groups, targets, d=d_b, bounds={"L": 0.0, "U": 1.5}, family="kl")
check("bound_viol_count", bv["bound_viol_count"], 2)
check("bound_viol_max", bv["bound_viol_max"], 0.5)
check("bound_viol_mean", bv["bound_viol_mean"], 1 / 6)

# -----------------------------------------------------------------------------
# Case (c): minimax objective-value agreement (Blocker G) -- weight-vector
# correlation must NOT be computed for the minimax family.
# -----------------------------------------------------------------------------
print("\n== Case (c): minimax objective-value agreement (Blocker G) ==")
w_mm = np.array([0.8, 1.1, 0.9, 1.3, 0.95, 0.95])
w_mm_ref = np.array([0.7, 1.2, 1.0, 1.2, 0.95, 0.95])  # deliberately different vector, same L-inf optimum face
ag = agreement(w_mm, w_mm_ref, family="minimax", obj_val=0.123456, obj_val_ref=0.123457, obj_tol=1e-5)
check_true("minimax agreement mode == 'objective_value' (not weight_vector)",
           ag["mode"] == "objective_value")
check("minimax abs_diff", ag["abs_diff"], 1.000000000001e-06, tol=1e-12)
check_true("minimax agree True within obj_tol", ag["agree"])
check_true("minimax agreement record has NO pearson/spearman fields (Blocker G)",
           "pearson" not in ag and "spearman" not in ag)

ag_disagree = agreement(w_mm, w_mm_ref, family="minimax", obj_val=0.05, obj_val_ref=0.20, obj_tol=1e-6)
check_true("minimax agree False when objective values differ beyond tol", not ag_disagree["agree"])

# w0 is uniform (sd=0) which would make corrcoef emit a NaN warning; use w_b.
ag_kl = agreement(w_b, np.array([1.1, 0.9, 1.0, 1.0, 2.05, 1.95]), family="kl")
check_true("kl-family agreement mode == 'weight_vector'", ag_kl["mode"] == "weight_vector")
check_true("kl-family agreement has pearson field", "pearson" in ag_kl)

# -----------------------------------------------------------------------------
print(f"\n{'GOLDEN PASS' if fail_count == 0 else 'GOLDEN FAIL'}: {fail_count} assertion(s) failed.")
if __name__ == "__main__":
    sys.exit(0 if fail_count == 0 else 1)
