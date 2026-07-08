"""benchmarks/study/common/metrics.py

Shared quality metrics for the leafblower benchmark study (DESIGN.md Section 6).
Python twin of metrics.R -- kept formula-identical for the R<->Python rtol=1e-6
parity gate. Reconciles THREE existing, mutually-inconsistent implementations:
  benchmarks/stepstone_all_methods.R:22   fit_metrics      (:34 t>0&S>0 gate BUG, :37 d_i=1 BUG)
  benchmarks/python_ipf_benchmark.py:31   compute_metrics  (same two bugs)
  benchmarks/allmethod_bench.R:30         compute_metrics  (3rd, differently-normalised weight_kl/chi2)
None of the three files above are modified here (read-only reference).

Guards enforced (see metrics.R docstring for full rationale):
 - margin KL is divergent-aware: T_kj>0 & p_kj<=0 -> +inf, `divergent=True`,
   term NEVER dropped (no t>0&p>0 gate).
 - weight_kl uses the REAL per-problem design weights d_i (never hardcoded 1),
   reported family-native with an explicit neutral-axis caveat.
 - DEFF reported as "Kish weighting DEFF/UWE" (raw w) plus a separate
   g-weight efficiency 1+CV^2(g), g=w/d, when d_i != 1 (Deville-Saerndal).
 - No cancellation: ratio-of-sums-of-squares forms only, never
   mean(x**2)-mean(x)**2. Zero-design-weight rows excluded from weight_kl/g.
 - RQ5 agreement: weight-vector stats only for strictly-convex families;
   minimax judged on achieved-L-infinity objective-value agreement
   (Blocker G -- non-unique L-inf LP optimum).

Single-thread BLAS is forced BEFORE numpy import for R<->Python parity
determinism (project convention; see /home/dd/Gemini/leafblower/CLAUDE.md).
"""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import math  # noqa: E402
import numpy as np  # noqa: E402

KL_STRICTLY_CONVEX_FAMILIES = {"kl", "chi2", "logit"}
KL_NATIVE_FAMILIES = {"kl", "raking", "sinkhorn", "greenkhorn", "oris", "oris_soft", "newton_kl"}


# ---- margin KL / L-inf / L1 -------------------------------------------------

def margin_one(w, group, target):
    """One margin's proportions, KL term, and error stats.

    w: 1D array-like, length n.
    group: 1D array-like, length n, category label per obs for this margin.
    target: dict {level: target proportion}, values sum to 1.
    """
    w = np.asarray(w, dtype=float)
    group = np.asarray(group)
    W = w.sum()
    levels = list(target.keys())
    S = np.array([w[group == lv].sum() for lv in levels], dtype=float)
    p = S / W
    T = np.array([target[lv] for lv in levels], dtype=float)
    err = np.abs(p - T)
    term = np.empty_like(T)
    for i in range(len(T)):
        if T[i] == 0:
            term[i] = 0.0
        elif p[i] <= 0:
            term[i] = math.inf
        else:
            term[i] = T[i] * math.log(T[i] / p[i])
    return dict(
        kl=float(term.sum()),
        divergent=bool(np.any((T > 0) & (p <= 0))),
        linf=float(err.max()), l1=float(err.sum()),
        p=p, T=T,
    )


def margin_stats(w, groups, targets):
    """Margin KL (mean/max across margins), divergence flag, margin L-inf/L1.

    groups: dict {margin_name: category array}, one entry per margin.
    targets: dict {margin_name: {level: target proportion}}, same keys as groups.
    """
    ks = list(targets.keys())
    per = {k: margin_one(w, groups[k], targets[k]) for k in ks}
    kl_k = np.array([per[k]["kl"] for k in ks])
    divergent = np.array([per[k]["divergent"] for k in ks])
    linf_k = np.array([per[k]["linf"] for k in ks])
    l1_k = np.array([per[k]["l1"] for k in ks])
    return dict(
        marg_kl_mean=float(kl_k.mean()), marg_kl_max=float(kl_k.max()),
        divergent=bool(divergent.any()),
        divergent_margins=[k for k in ks if per[k]["divergent"]],
        margin_linf=float(linf_k.max()), margin_l1=float(l1_k.sum()),
        per_margin=per,
    )


# ---- ESS / Kish DEFF / g-weight efficiency ---------------------------------

def _deff_ratio(x):
    """Ratio-of-sums-of-squares form of 1+CV^2(x) == n*sum(x^2)/sum(x)^2.

    Cancellation-free: no mean(x^2)-mean(x)^2 subtraction of close numbers.
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    return n * np.sum(x ** 2) / np.sum(x) ** 2


def ess_deff(w, d=None):
    """ESS (Kish), Kish weighting DEFF/UWE (raw w), g-weight efficiency
    1+CV^2(g) (g=w/d) when design weights are non-uniform.
    """
    w = np.asarray(w, dtype=float)
    n = len(w)
    ess = float(w.sum() ** 2 / np.sum(w ** 2))
    out = dict(
        ESS=ess, n=n, DEFF_kish=n / ess,
        DEFF_kish_label="Kish weighting DEFF/UWE (raw w) -- NOT the true design effect",
        g_weighted=False, DEFF_g=float("nan"),
    )
    if d is not None:
        d = np.asarray(d, dtype=float)
        use = d > 0
        if np.any(~use):
            out["excluded_zero_design"] = int(np.sum(~use))
        if np.any(d[use] != 1):
            g = w[use] / d[use]
            out["DEFF_g"] = float(_deff_ratio(g))
            out["g_weighted"] = True
    return out


# ---- closeness-to-design weight_kl (family-native) -------------------------

def weight_kl(w, d, family=None):
    """weight_kl = sum(w_i*log(w_i/d_i)), REAL per-problem d_i (never hardcoded 1).

    Rows with d_i==0 excluded (0*inf guard on a divergent shift); w_i==0
    contributes 0 by convention (0*log(0/d)=0).
    """
    w = np.asarray(w, dtype=float)
    d = np.asarray(d, dtype=float)
    use = d > 0
    w_u, d_u = w[use], d[use]
    safe_w = np.where(w_u > 0, w_u, 1.0)  # avoid log(0); masked out by the where() below
    term = np.where(w_u > 0, w_u * np.log(safe_w / d_u), 0.0)
    native = family in KL_NATIVE_FAMILIES
    note = ("weight_kl IS this family's objective (KL/raking dual) -- not a neutral cross-family axis"
            if native else
            "weight_kl reported as closeness-to-design diagnostic (non-native for this family)")
    return dict(weight_kl=float(term.sum()), excluded_zero_design=int(np.sum(~use)),
                family=family, neutral_axis_note=note)


# ---- bound violation --------------------------------------------------------

def bound_violation(w, L=0.0, U=math.inf):
    """Count and max/mean magnitude of max(0, L-w, w-U). L/U scalar or per-obs array."""
    w = np.asarray(w, dtype=float)
    viol = np.maximum(0.0, np.maximum(L - w, w - U))
    return dict(count=int(np.sum(viol > 0)), max=float(viol.max()), mean=float(viol.mean()))


# ---- RQ5 agreement (Blocker G: minimax excluded from vector correlation) --

def _rank(x):
    """Average-tie ranks (1-based), matching R's cor(method="spearman")."""
    x = np.asarray(x, dtype=float)
    order = np.argsort(x, kind="mergesort")
    ranks = np.empty(len(x), dtype=float)
    ranks[order] = np.arange(1, len(x) + 1)
    sorted_x = x[order]
    i = 0
    while i < len(x):
        j = i
        while j + 1 < len(x) and sorted_x[j + 1] == sorted_x[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = ranks[order[i:j + 1]].mean()
        i = j + 1
    return ranks


def agreement(w, w_ref, family, obj_val=None, obj_val_ref=None, obj_tol=1e-6):
    """Weight-vector agreement for strictly-convex families; achieved-L-inf
    objective-value agreement for minimax (unique-optimum caveat, Blocker G).
    """
    if family in KL_STRICTLY_CONVEX_FAMILIES:
        w = np.asarray(w, dtype=float)
        w_ref = np.asarray(w_ref, dtype=float)
        pearson = float(np.corrcoef(w, w_ref)[0, 1])
        spearman = float(np.corrcoef(_rank(w), _rank(w_ref))[0, 1])
        cosine = float(np.dot(w, w_ref) / (np.linalg.norm(w) * np.linalg.norm(w_ref)))
        return dict(family=family, mode="weight_vector", pearson=pearson, spearman=spearman,
                    max_abs_diff=float(np.max(np.abs(w - w_ref))), cosine=cosine)
    if obj_val is None or obj_val_ref is None:
        raise ValueError("minimax/non-strictly-convex agreement requires obj_val and obj_val_ref")
    d = abs(obj_val - obj_val_ref)
    return dict(family=family, mode="objective_value", obj_val=obj_val, obj_val_ref=obj_val_ref,
                abs_diff=d, agree=bool(d <= obj_tol))


# ---- top-level convenience wrapper -----------------------------------------

def compute_metrics(w, groups, targets, d=None, bounds=None, family=None,
                     w_ref=None, obj_val=None, obj_val_ref=None, obj_tol=1e-6):
    """Full metrics record for one (weights, problem) pair. See metrics.R for
    the parameter contract (kept identical across languages)."""
    w = np.asarray(w, dtype=float)
    if d is None:
        d = np.ones_like(w)
    ms = margin_stats(w, groups, targets)
    ed = ess_deff(w, d)
    wk = weight_kl(w, d, family)
    out = dict(
        n=len(w), W=float(w.sum()),
        marg_kl_mean=ms["marg_kl_mean"], marg_kl_max=ms["marg_kl_max"],
        marg_kl_divergent=ms["divergent"], marg_kl_divergent_margins=ms["divergent_margins"],
        margin_linf=ms["margin_linf"], margin_l1=ms["margin_l1"],
        ESS=ed["ESS"], DEFF_kish=ed["DEFF_kish"], DEFF_kish_label=ed["DEFF_kish_label"],
        DEFF_g=ed["DEFF_g"], g_weighted=ed["g_weighted"],
        weight_kl=wk["weight_kl"], weight_kl_excluded_zero_design=wk["excluded_zero_design"],
        weight_kl_note=wk["neutral_axis_note"],
        wmin=float(w.min()), wmed=float(np.median(w)), wmax=float(w.max()),
    )
    if bounds is not None:
        bv = bound_violation(w, bounds.get("L", 0.0), bounds.get("U", math.inf))
        out["bound_viol_count"] = bv["count"]
        out["bound_viol_max"] = bv["max"]
        out["bound_viol_mean"] = bv["mean"]
    if w_ref is not None:
        out["agreement"] = agreement(w, w_ref, family, obj_val, obj_val_ref, obj_tol)
    return out
