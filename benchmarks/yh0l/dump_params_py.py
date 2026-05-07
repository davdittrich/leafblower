#!/usr/bin/env python3
"""T1: dump all params dict entries for ieppa_soft as Python resolves them.

Replicates _harvest.py param resolution for:
  harvest(method='ieppa_soft', convergence={'tol': 1e-4},
          max_weight=5, max_iterations=3000)
Does NOT call the solver.

Also documents C++ binding default for each key NOT in the dict (from _bindings.cpp).
"""

import sys
import os

# Use the installed leafblower package (main repo install; worktree shares .so).
# If worktree python/ has a built .so, prefer it; else fall back to installed.
_wt = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_wt_py = os.path.join(_wt, "python")
if os.path.exists(os.path.join(_wt_py, "leafblower", "_leafblower.cpython-314-x86_64-linux-gnu.so")) or \
   os.path.exists(os.path.join(_wt_py, "leafblower", "_leafblower.cpython-312-x86_64-linux-gnu.so")):
    sys.path.insert(0, _wt_py)
# else: use installed package

from leafblower._harvest import (
    _parse_convergence, _parse_sor, _resolve_sor,
    _METRIC_MAP,
)

# ---- Call parameters matching the test invocation ----
method         = "ieppa_soft"
convergence    = {"tol": 1e-4}
min_weight     = 0.0
max_weight     = 5.0
max_iterations = 3000
verbose        = 0
capacity_penalty = None
alm_penalty    = None
bounds_mode    = "cell"
homotopy_levels       = 1
homotopy_start_factor = 1.0
homotopy_end_factor   = 1.0
homotopy_budget_p     = 0.5
scheduler      = "round_robin"
eta_schedule   = "fixed"
eta_start      = 1.0
eta_end        = 1.0
eta_schedule_power = 0.5
accelerate     = False
newton_tsvd_ratio = 1e-8
alm_penalty    = None
ridge_lambda   = 0.0
sor            = None   # Python default: sor=None → _resolve_sor(None) → {} → _parse_sor({}) → disabled

# ---- Replicate resolution logic ----

# convergence = {"tol": 1e-4}:
#   No shorthand key → default_metric="max_err", default_rule="improvement", default_tol=0.001
#   tol overridden to 1e-4 by conv.get("tol", default_tol)
#   rule="improvement" → pct_tol=tol=1e-4, absolute_tol=0.0
pct_tol, absolute_tol, metric, rule, stop_when = _parse_convergence(convergence)

# Per-method metric override (mirrors R harvest.R:308-331)
# No explicit metric in convergence dict → ieppa_soft → marginal_kl (int 6)
_conv = convergence or {}
if not any(k in _conv for k in ("metric", "improvement", "pct", "absolute")):
    _method_metric_map = {
        "ieppa": "marginal_kl", "ieppa_soft": "marginal_kl",
        "raking": "kl", "greenkhorn": "kl",
        "sinkhorn": "kl", "newton_kl": "kl",
        "greg": "chi2",
    }
    _m = _method_metric_map.get(method.lower())
    if _m is not None:
        metric = _METRIC_MAP[_m]

# SOR
sor_enabled, sor_auto, sor_omega_init, sor_omega_min, sor_omega_fixed, sor_burnin = _parse_sor(
    _resolve_sor(sor)
)

# algorithm int
alg_map = {
    "ieppa": 1, "raking": 3,
    "sinkhorn": 4, "chebyshev": 5, "greg": 6,
    "ieppa_soft": 8, "greenkhorn": 9, "logit": 10, "newton_kl": 11,
}
alg_int = alg_map[method.lower()]

# bounds_mode int
_bounds_mode_int = {"cell": 0, "unit": 1}[bounds_mode]

# scheduler int
_scheduler_int = {"round_robin": 0, "greedy": 1}[scheduler]

# eta_schedule int
_eta_mode_int = {"fixed": 0, "tang_dynamic": 1}[eta_schedule]

# tol_abs (legacy): rule="improvement" → absolute_tol=0.0, use 1e-6 fallback
tol_abs_legacy = absolute_tol if absolute_tol > 0.0 else 1e-6

# Build params dict (exactly as _harvest.py does)
params = {
    "min_weight":     min_weight,
    "max_weight":     max_weight,
    "inner_max_iter": max_iterations,
    "outer_max_iter": max_iterations,   # mirrors R bridge: outer = inner = max_iterations
    "tol_abs":        tol_abs_legacy,
    "verbose":        verbose,
    "algorithm":      alg_int,
    "epsilon":        0.05,
    "bounds_mode":    _bounds_mode_int,
    # Convergence config (WU-G)
    "pct_tol":        pct_tol,
    "absolute_tol":   absolute_tol,
    "metric":         metric,
    "rule":           rule,
    "stop_when":      stop_when,
    # SOR config (WU-F)
    "sor_enabled":    sor_enabled,
    "sor_auto":       sor_auto,
    "sor_omega_init": sor_omega_init,
    "sor_omega_min":  sor_omega_min,
    "sor_omega_fixed": sor_omega_fixed,
    "sor_burnin":     sor_burnin,
    # Homotopy config (PY-1)
    "homotopy_levels":       homotopy_levels,
    "homotopy_start_factor": homotopy_start_factor,
    "homotopy_end_factor":   homotopy_end_factor,
    "homotopy_budget_p":     homotopy_budget_p,
    # Scheduler / eta (PY-1)
    "scheduler":             _scheduler_int,
    "eta_mode":              _eta_mode_int,
    "eta_start":             eta_start,
    "eta_end":               eta_end,
    "eta_schedule_power":    eta_schedule_power,
    # Method-specific (PY-1)
    "newton_tsvd_ratio":     newton_tsvd_ratio,
    # SRAA / ALM (PY-2)
    "accelerate":            int(accelerate),
    "ridge_lambda":          ridge_lambda,
}
# capacity_penalty and alm_penalty are conditionally added
if capacity_penalty is not None:
    params["capacity_penalty"] = capacity_penalty
if alm_penalty is not None and alm_penalty > 0.0:
    params["alm_penalty"] = alm_penalty

# C++ binding defaults (from _bindings.cpp + rk_params_init in c_api.cpp).
# For keys NOT in params dict, the rk_params_init default fires.
# These are verified from: src/c_api.cpp:rk_params_init (lines 51-89)
CPP_DEFAULTS = {
    # Keys present in params dict (no default fires)
    "min_weight":            ("IN_DICT",  min_weight),
    "max_weight":            ("IN_DICT",  max_weight),
    "inner_max_iter":        ("IN_DICT",  max_iterations),
    "outer_max_iter":        ("IN_DICT",  max_iterations),
    "tol_abs":               ("IN_DICT",  tol_abs_legacy),
    "verbose":               ("IN_DICT",  verbose),
    "algorithm":             ("IN_DICT",  alg_int),
    "epsilon":               ("IN_DICT",  0.05),
    "bounds_mode":           ("IN_DICT",  _bounds_mode_int),
    "pct_tol":               ("IN_DICT",  pct_tol),
    "absolute_tol":          ("IN_DICT",  absolute_tol),
    "metric":                ("IN_DICT",  metric),
    "rule":                  ("IN_DICT",  rule),
    "stop_when":             ("IN_DICT",  stop_when),
    "sor_enabled":           ("IN_DICT",  sor_enabled),
    "sor_auto":              ("IN_DICT",  sor_auto),
    "sor_omega_init":        ("IN_DICT",  sor_omega_init),
    "sor_omega_min":         ("IN_DICT",  sor_omega_min),
    "sor_omega_fixed":       ("IN_DICT",  sor_omega_fixed),
    "sor_burnin":            ("IN_DICT",  sor_burnin),
    "homotopy_levels":       ("IN_DICT",  homotopy_levels),
    "homotopy_start_factor": ("IN_DICT",  homotopy_start_factor),
    "homotopy_end_factor":   ("IN_DICT",  homotopy_end_factor),
    "homotopy_budget_p":     ("IN_DICT",  homotopy_budget_p),
    "scheduler":             ("IN_DICT",  _scheduler_int),
    "eta_mode":              ("IN_DICT",  _eta_mode_int),
    "eta_start":             ("IN_DICT",  eta_start),
    "eta_end":               ("IN_DICT",  eta_end),
    "eta_schedule_power":    ("IN_DICT",  eta_schedule_power),
    "newton_tsvd_ratio":     ("IN_DICT",  newton_tsvd_ratio),
    "accelerate":            ("IN_DICT",  int(accelerate)),
    "ridge_lambda":          ("IN_DICT",  ridge_lambda),
    # Keys NOT in params dict (capacity_penalty=None, alm_penalty=None):
    # C++ binding uses if(params_dict.contains(key)) pattern.
    # When absent, rk_params_init leaves field at its default.
    # capacity_penalty: rk_params_init does memset(p,0) → 0.0
    #   (field p->capacity_penalty is NOT set in rk_params_init; memset=0 → 0.0)
    #   src/c_api.cpp line 51: memset(p, 0, ...) covers this field.
    #   Effective C++ default: 0.0 (treated as <=0.0 → auto = M_cell/n in ieppa solver)
    "capacity_penalty":      ("OMITTED_CPP_DEFAULT", 0.0),  # c_api.cpp:memset→0.0; ieppa_soft treats <=0 as auto
    # alm_penalty: same — not set in rk_params_init, memset=0 → 0.0
    #   src/c_api.cpp line 51: memset covers this.
    #   Effective C++ default: 0.0 (0.0 = inactive per struct comment)
    "alm_penalty":           ("OMITTED_CPP_DEFAULT", 0.0),  # c_api.cpp:memset→0.0; 0.0=inactive
}

print("=== Python param dump: harvest(method='ieppa_soft', convergence={'tol': 1e-4}, max_weight=5, max_iterations=3000) ===\n")

# Slots 1-3 and slot 10: data-dependent args — structural summary (not in params dict; passed directly to C++)
print("--- data/input args (positional, data-dependent; structural summary) ---")
print("  group_ids_r             = VECSXP[K]; each INTSXP len=n (0-idx codes, -1=NA). K=n_margins, n=nrow(data)")
print("  cat_counts_r            = INTSXP[K]; per-margin level counts. sum=total_cells")
print("  targets_r               = VECSXP[K]; each REALSXP len=cat_counts[k]. sum(targets[k])≈1")
print("  n_obs                   = INTSXP scalar = nrow(data)")
print("  sw_vec (start_weights)  = REALSXP len=n or NULL → uniform 1.0")

print("\n--- params dict (passed to C++ calibrate()) ---")
for k, v in params.items():
    status = "IN_DICT"
    cpp_default = CPP_DEFAULTS.get(k, ("IN_DICT", None))[1]
    print(f"  {k:<26} = {repr(v)}")

print("\n--- Keys OMITTED from params dict (capacity_penalty, alm_penalty) ---")
print(f"  {'capacity_penalty':<26}: NOT in dict (capacity_penalty=None)")
print(f"    → _bindings.cpp:131 if(params_dict.contains('capacity_penalty')) NOT taken")
print(f"    → C++ field p.capacity_penalty left at rk_params_init default: 0.0 (memset, c_api.cpp:51)")
print(f"    → ieppa_soft treats <=0.0 as auto (M_cell/n from cell_table)")
print(f"  {'alm_penalty':<26}: NOT in dict (alm_penalty=None)")
print(f"    → _bindings.cpp:136 if(params_dict.contains('alm_penalty')) NOT taken")
print(f"    → C++ field p.alm_penalty left at rk_params_init default: 0.0 (memset, c_api.cpp:51)")
print(f"    → 0.0 = inactive (struct comment leafblower.h:100)")

print("\n--- Convergence resolution detail ---")
print(f"  convergence input      = {convergence!r}")
print(f"  pct_tol (resolved)     = {pct_tol!r}  (tol=1e-4 with rule=improvement)")
print(f"  absolute_tol           = {absolute_tol!r}")
print(f"  metric (int)           = {metric}  (marginal_kl=6, per ieppa_soft method override)")
print(f"  rule (int)             = {rule}  (improvement=1)")
print(f"  stop_when (int)        = {stop_when}  (any=0)")

print("\n--- SOR resolution ---")
print(f"  sor=None → _resolve_sor(None) → {{}} → _parse_sor({{}}): DISABLED")
print(f"  sor_enabled  = {sor_enabled}")
print(f"  sor_auto     = {sor_auto}")
print(f"  sor_omega_init={sor_omega_init}")
print(f"  sor_omega_min ={sor_omega_min}")
print(f"  sor_omega_fixed={sor_omega_fixed}")
print(f"  sor_burnin   = {sor_burnin}")

print("\nDone.")
