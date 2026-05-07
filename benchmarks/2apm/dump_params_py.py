#!/usr/bin/env python3
"""T1: dump all params dict entries for chebyshev as Python resolves them.

Replicates _harvest.py param resolution for:
  harvest(method='chebyshev', convergence={'tol': 1e-4},
          max_weight=5, max_iterations=3000)
Does NOT call the solver.

Also documents C++ binding default for each key NOT in the dict (from _bindings.cpp).
"""

import sys
import os

if not os.path.exists("benchmarks/2apm"):
    sys.exit(f"Run from project root. Current wd={os.getcwd()}")

# Prefer worktree python/ build if present; else fall back to installed package.
# Resolve from cwd (project root), not __file__, so worktree and main repo both work.
_wt_py = os.path.join(os.getcwd(), "python")
if os.path.isdir(_wt_py) and (
    os.path.exists(os.path.join(_wt_py, "leafblower", "_leafblower.cpython-314-x86_64-linux-gnu.so")) or
    os.path.exists(os.path.join(_wt_py, "leafblower", "_leafblower.cpython-312-x86_64-linux-gnu.so"))
):
    sys.path.insert(0, _wt_py)

from leafblower._harvest import (
    _parse_convergence, _parse_sor, _resolve_sor,
    _METRIC_MAP,
)

method         = "chebyshev"
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
ridge_lambda   = 0.0
sor            = None   # Python default: sor=None → _resolve_sor(None) → {} → _parse_sor({}) → disabled

# ---- Replicate resolution logic ----

# convergence = {"tol": 1e-4}:
# No shorthand key → rule=improvement, pct_tol=tol=1e-4, absolute_tol=0.0
pct_tol, absolute_tol, metric, rule, stop_when = _parse_convergence(convergence)

# Per-method metric override (mirrors R harvest.R:308-331)
# chebyshev is NOT in the override map → no override → metric stays at parse_convergence default
# (Python comment: "chebyshev: omitted — max_err is its natural objective (L-inf)")
# _parse_convergence default for no explicit metric = max_err (int 0)
_conv = convergence or {}
if not any(k in _conv for k in ("metric", "improvement", "pct", "absolute")):
    _method_metric_map = {
        "ieppa": "marginal_kl", "ieppa_soft": "marginal_kl",
        "raking": "kl", "greenkhorn": "kl",
        "sinkhorn": "kl", "newton_kl": "kl",
        "greg": "chi2",
        # chebyshev: omitted — max_err is its natural objective (L-inf)
    }
    _m = _method_metric_map.get(method.lower())
    if _m is not None:
        metric = _METRIC_MAP[_m]
# chebyshev not in map → metric stays at default from _parse_convergence = 0 (max_err)

sor_enabled, sor_auto, sor_omega_init, sor_omega_min, sor_omega_fixed, sor_burnin = _parse_sor(
    _resolve_sor(sor)
)

alg_map = {
    "ieppa": 1, "raking": 3,
    "sinkhorn": 4, "chebyshev": 5, "greg": 6,
    "ieppa_soft": 8, "greenkhorn": 9, "logit": 10, "newton_kl": 11,
}
alg_int = alg_map[method.lower()]  # 5

_bounds_mode_int = {"cell": 0, "unit": 1}[bounds_mode]
_scheduler_int = {"round_robin": 0, "greedy": 1}[scheduler]
_eta_mode_int = {"fixed": 0, "tang_dynamic": 1}[eta_schedule]

# tol_abs (legacy): rule="improvement" → absolute_tol=0.0, use 1e-6 fallback
tol_abs_legacy = absolute_tol if absolute_tol > 0.0 else 1e-6

# ---- Build params dict (exactly as _harvest.py does) ----
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
    "pct_tol":        pct_tol,
    "absolute_tol":   absolute_tol,
    "metric":         metric,
    "rule":           rule,
    "stop_when":      stop_when,
    "sor_enabled":    sor_enabled,
    "sor_auto":       sor_auto,
    "sor_omega_init": sor_omega_init,
    "sor_omega_min":  sor_omega_min,
    "sor_omega_fixed": sor_omega_fixed,
    "sor_burnin":     sor_burnin,
    "homotopy_levels":       homotopy_levels,
    "homotopy_start_factor": homotopy_start_factor,
    "homotopy_end_factor":   homotopy_end_factor,
    "homotopy_budget_p":     homotopy_budget_p,
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
# capacity_penalty OMITTED from dict → C++ binding skips; field left at memset=0.0
# chebyshev: alm block gated by st.use_admm_capacity; capacity_mu harmless for non-ieppa_soft
if capacity_penalty is not None:
    params["capacity_penalty"] = capacity_penalty
# alm_penalty OMITTED from dict → C++ binding skips; field left at memset=0.0 (inactive)
if alm_penalty is not None and alm_penalty > 0.0:
    params["alm_penalty"] = alm_penalty

# ---- C++ binding defaults (from _bindings.cpp + rk_params_init in c_api.cpp). ----
# For keys NOT in params dict, the rk_params_init default fires.
# These are verified from: src/c_api.cpp:rk_params_init (lines 51-89)
CPP_DEFAULTS = {
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
    "capacity_penalty":      ("OMITTED_CPP_DEFAULT", 0.0),  # c_api.cpp:memset→0.0; chebyshev ALM gated
    "alm_penalty":           ("OMITTED_CPP_DEFAULT", 0.0),  # c_api.cpp:memset→0.0; 0.0=inactive
}

print("=== Python param dump: harvest(method='chebyshev', convergence={'tol': 1e-4}, max_weight=5, max_iterations=3000) ===\n")
print("--- data/input args (positional, data-dependent; structural summary) ---")
print("  group_ids_r             = VECSXP[K]; each INTSXP len=n (0-idx codes, -1=NA). K=n_margins, n=nrow(data)")
print("  cat_counts_r            = INTSXP[K]; per-margin level counts. sum=total_cells")
print("  targets_r               = VECSXP[K]; each REALSXP len=cat_counts[k]. sum(targets[k])≈1")
print("  n_obs                   = INTSXP scalar = nrow(data)")
print("  sw_vec (start_weights)  = REALSXP len=n or NULL → uniform 1.0")

print("\n--- params dict (passed to C++ calibrate()) ---")
for k, v in params.items():
    print(f"  {k:<26} = {repr(v)}")

print("\n--- Keys OMITTED from params dict (capacity_penalty, alm_penalty) ---")
print("OMITTED_KEY: capacity_penalty | C_DEFAULT: 0.0 | CITATION: c_api.cpp:51")
print("OMITTED_KEY: alm_penalty | C_DEFAULT: 0.0 | CITATION: c_api.cpp:51")

print("\n--- Convergence resolution detail ---")
print(f"  convergence input      = {convergence!r}")
print(f"  pct_tol (resolved)     = {pct_tol!r}  (tol=1e-4 with rule=improvement)")
print(f"  absolute_tol           = {absolute_tol!r}")
print(f"  metric (int)           = {metric}  (max_err=0; chebyshev has no override — L-inf natural objective)")
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
