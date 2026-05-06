from __future__ import annotations
import math
import warnings
import numpy as np
from typing import TYPE_CHECKING, Dict, Optional

if TYPE_CHECKING:
    import pandas as pd  # only for type checker; runtime import below

try:
    import pandas as pd  # type: ignore[no-redef]
    _PANDAS_AVAILABLE = True
except ImportError:
    _PANDAS_AVAILABLE = False

from ._leafblower import calibrate

_TARGET_SUM_TOL = 5e-7  # matches R harvest() tolerance for target proportions summing to 1

_METRIC_MAP = {
    "max_err": 0, "mean_err": 1, "kl": 2, "chi2": 3,
    "grake_norm": 4, "l1_weight": 5, "marginal_kl": 6,
}
_RULE_MAP = {"threshold": 0, "improvement": 1, "plateau": 2}
_STOP_WHEN_MAP = {"any": 0, "all": 1}
_KNOWN_CONVERGENCE_KEYS = frozenset({
    "metric", "rule", "tol", "pct", "absolute", "improvement", "stop_when",
})
# Reverse-lookup arrays mirror CalibMetric/CalibRule enum order in leafblower.h
_METRIC_NAMES = ["max_err", "mean_err", "kl", "chi2", "grake_norm", "l1_weight", "marginal_kl"]
_RULE_NAMES   = ["threshold", "improvement", "plateau"]


def _parse_convergence(conv):
    """Derive pct_tol, absolute_tol, metric, rule, stop_when from convergence dict.

    Mirrors R's parse_convergence() logic exactly:
      - Default: max_err + improvement + tol=0.001
      - 'improvement' key: max_err + improvement rule, tol = improvement value
      - 'pct' key: l1_weight + plateau rule, tol = pct value
      - 'absolute' key alone: max_err + threshold rule, tol = absolute value
      - 'metric'/'rule'/'tol' keys override the derived defaults
    """
    if conv is None:
        conv = {}
    unknown = set(conv) - _KNOWN_CONVERGENCE_KEYS
    if unknown:
        raise ValueError(f"unknown convergence key(s): {', '.join(sorted(unknown))}")

    explicit_impr = "improvement" in conv
    explicit_pct  = "pct" in conv
    explicit_abs  = "absolute" in conv

    if explicit_impr and explicit_abs and "stop_when" not in conv:
        raise ValueError(
            "convergence: 'improvement' and 'absolute' cannot be combined without "
            "'stop_when'. Use stop_when='any' (fire on either) or 'all' (require both).")

    if explicit_impr:
        default_metric, default_rule, default_tol = "max_err", "improvement", float(conv["improvement"])
    elif explicit_pct:
        default_metric, default_rule, default_tol = "l1_weight", "plateau", float(conv["pct"])
    elif not explicit_abs:
        default_metric, default_rule, default_tol = "max_err", "improvement", 0.001
    else:
        default_metric, default_rule, default_tol = "max_err", "threshold", float(conv["absolute"])

    metric_str    = conv.get("metric", default_metric)
    rule_str      = conv.get("rule", default_rule)
    tol           = float(conv.get("tol", default_tol))
    abs_tol_raw   = float(conv.get("absolute", 0.0))
    stop_when_str = conv.get("stop_when", "any")

    if metric_str not in _METRIC_MAP:
        raise ValueError(f"metric must be one of {list(_METRIC_MAP)}")
    if rule_str not in _RULE_MAP:
        raise ValueError(f"rule must be one of {list(_RULE_MAP)}")
    if stop_when_str not in _STOP_WHEN_MAP:
        raise ValueError(f"stop_when must be 'any' or 'all'")

    if rule_str == "threshold":
        pct_tol      = 0.0
        absolute_tol = tol
    else:
        pct_tol      = tol
        absolute_tol = abs_tol_raw

    return pct_tol, absolute_tol, _METRIC_MAP[metric_str], _RULE_MAP[rule_str], _STOP_WHEN_MAP[stop_when_str]


def _resolve_sor(sor):
    """Map sor=None to R default (auto-enabled). sor=False/sor={} = explicit disable."""
    if sor is None:
        return {"auto": True, "omega_min": 0.3}  # R default: list(auto=TRUE, omega_min=0.3)
    if sor is False:
        return {}  # explicit disable
    return sor  # passthrough dict


def _parse_sor(sor):
    """Mirror R parse_sor(): returns (enabled, auto, omega_init, omega_min, omega_fixed, burnin)."""
    if not sor:  # empty dict = disabled (from _resolve_sor(False))
        return 0, 0, 1.0, 0.3, -1.0, 20
    enabled = 1
    auto = 1 if sor.get("auto", True) else 0
    omega_init = float(sor.get("omega_init", 1.0))
    omega_min = float(sor.get("omega_min", 0.3))
    omega_fixed = float(sor.get("omega", -1.0))
    burnin = int(sor.get("burnin", 20))
    return enabled, auto, omega_init, omega_min, omega_fixed, burnin


def harvest(
    data,
    targets: Dict[str, Dict[str, float]],
    min_weight: float = 0.0,
    max_weight: float = 5.0,
    method: str = "ieppa",
    verbose: int = 0,
    max_iterations: int = 500,
    start_weights: Optional[np.ndarray] = None,
    attach_weights: bool = True,
    weight_column: str = "weights",
    convergence: Optional[Dict] = None,
    sor: Optional[Dict] = None,
    bounds_mode: str = "cell",
    homotopy_levels: int = 1,
    homotopy_start_factor: float = 1.0,
    homotopy_end_factor: float = 1.0,
    homotopy_budget_p: float = 0.5,
    scheduler: str = "round_robin",
    eta_schedule: str = "fixed",
    eta_start: float = 1.0,
    eta_end: float = 1.0,
    eta_schedule_power: float = 0.5,
    capacity_penalty: Optional[float] = None,
    newton_tsvd_ratio: float = 1e-8,
    accelerate: bool = False,
    alm_penalty: Optional[float] = None,
    add_na_proportion: bool = False,
    auto_collapse: bool = False,
    collapse_vars=None,
    design_weights=None,
    **_kwargs,  # absorbed for forward-compat; not passed to R
):
    """
    Calibrate survey weights. Drop-in for R leafblower::harvest().

    Parameters
    ----------
    data : pd.DataFrame or dict of lists
    targets : dict of dicts, e.g. {"age": {"18-34": 0.3, "35+": 0.7}}
    min_weight : float, lower bound on weights (default 0 = no bound)
    max_weight : float, upper bound on weights (default 5)
    method : str, one of "ieppa" (default), "raking", "sinkhorn",
        "chebyshev", "greg", "ieppa_soft", "greenkhorn", "logit", "newton_kl"
    verbose : int, 0=silent, 1=progress, 2=debug
    max_iterations : int, inner BCD max sweeps per outer iter (default 500)
    start_weights : optional 1D float64 array of initial weights
    attach_weights : if True, return DataFrame with weights column appended
    weight_column : name of the weights column (default "weights")
    convergence : dict controlling the stopping criterion. Keys:
        "metric" (str) — one of "max_err" (default), "mean_err", "kl", "chi2",
            "grake_norm", "l1_weight", "marginal_kl".
        "rule" (str) — one of "improvement" (default), "threshold", "plateau".
        "tol" (float) — tolerance value (default 0.001).
        "pct" (float) — shorthand: activates l1_weight + plateau rule.
        "absolute" (float) — shorthand: activates max_err + threshold rule.
        "improvement" (float) — shorthand: activates max_err + improvement rule.
        "stop_when" (str) — "any" (default) or "all".
        Unknown keys raise ValueError.
    sor : dict for SOR adaptive under-relaxation (iEPPA only). None disables SOR. Keys:
        "auto" (bool, default True), "omega_min" (float, default 0.3),
        "omega" (float) for fixed omega, "burnin" (int, default 20).
    bounds_mode : str, "cell" (default) or "unit". "cell": per-cell aggregate bounds —
                  individual weights may fall outside [min_weight, max_weight] when base
                  weights are skewed within a cell. "unit": per-observation strict bounds
                  via intra-cell water-filling redistribution.
    Returns
    -------
    pd.DataFrame (if attach_weights=True) or np.ndarray
    """
    # Not-in-v1 hard stops (mirrors R harvest.R lines 244-249)
    if add_na_proportion is not False:
        raise ValueError("add_na_proportion is not supported in leafblower v1.")
    if auto_collapse is True:
        raise ValueError("auto_collapse is not supported in leafblower v1.")
    if collapse_vars is not None:
        raise ValueError("collapse_vars is not supported in leafblower v1.")

    # design_weights: alias for start_weights (mirrors R harvest.R lines 332-338)
    if design_weights is not None:
        if start_weights is not None:
            warnings.warn(
                "leafblower: both design_weights and start_weights supplied; design_weights ignored",
                UserWarning, stacklevel=2)
        else:
            start_weights = design_weights

    # Parse convergence and SOR
    pct_tol, absolute_tol, metric, rule, stop_when = _parse_convergence(convergence)
    sor_enabled, sor_auto, sor_omega_init, sor_omega_min, sor_omega_fixed, sor_burnin = _parse_sor(
        _resolve_sor(sor)
    )

    # Convert dict data to DataFrame; validate input type.
    if isinstance(data, dict):
        if not _PANDAS_AVAILABLE:
            raise ImportError("pandas required to use dict input; install with pip install pandas")
        data = pd.DataFrame(data)
    elif _PANDAS_AVAILABLE:
        if not isinstance(data, pd.DataFrame):
            raise TypeError("data must be a pd.DataFrame or dict")
    else:
        raise TypeError("data must be a dict when pandas is not installed")

    n = len(data)

    # Method mapping
    method_lc = method.lower()

    alg_map = {
        "ieppa": 1, "raking": 3,
        "sinkhorn": 4, "chebyshev": 5, "greg": 6,
        "ieppa_soft": 8, "greenkhorn": 9, "logit": 10, "newton_kl": 11,
    }  # "auto" (0) removed from Python user API; "grake" (7) removed; "lbfgsb" (2) removed — enum gap
    if method_lc not in alg_map:
        raise ValueError(f"method must be one of {list(alg_map)}")
    alg_int = alg_map[method_lc]

    if bounds_mode not in ("cell", "unit"):
        raise ValueError(f"bounds_mode must be 'cell' or 'unit', got {bounds_mode!r}")
    _bounds_mode_int = {"cell": 0, "unit": 1}[bounds_mode]

    # scheduler: str → int
    _scheduler_map = {"round_robin": 0, "greedy": 1}
    if scheduler not in _scheduler_map:
        raise ValueError(f"scheduler must be 'round_robin' or 'greedy', got {scheduler!r}")
    _scheduler_int = _scheduler_map[scheduler]

    # eta_schedule: str → int (maps to C field eta_mode)
    _eta_map = {"fixed": 0, "tang_dynamic": 1}
    if eta_schedule not in _eta_map:
        raise ValueError(f"eta_schedule must be 'fixed' or 'tang_dynamic', got {eta_schedule!r}")
    _eta_mode_int = _eta_map[eta_schedule]

    # Build group_ids and validate targets
    K = len(targets)
    group_ids_list = []
    cat_counts_list = []
    targets_list = []
    var_names = list(targets.keys())

    for varname in var_names:
        tgt_dict = targets[varname]
        levels = list(tgt_dict.keys())
        props = list(tgt_dict.values())
        if abs(sum(props) - 1.0) > _TARGET_SUM_TOL:
            raise ValueError(f"targets['{varname}'] does not sum to 1 (sum={sum(props):.8f})")

        col = data[varname]
        ncat = len(levels)

        if _PANDAS_AVAILABLE:
            # Vectorized encoding via pandas Categorical: O(n) in C, ~6x faster than
            # a Python for-loop. col.astype(str) handles mixed-type columns; pd.isna
            # entries map to codes=-1 automatically when observed=False.
            cat = pd.Categorical(col.astype(str).where(~pd.isna(col), other=np.nan),
                                 categories=[str(lv) for lv in levels])
            gid = cat.codes.astype(np.int32)  # -1 for NA/unknown levels
        else:
            level_to_idx = {lv: j for j, lv in enumerate(levels)}
            gid = np.empty(n, dtype=np.int32)
            for i, val in enumerate(col):
                if val is None or (isinstance(val, float) and np.isnan(val)):
                    gid[i] = -1
                else:
                    gid[i] = level_to_idx.get(str(val), -1)

        if len(gid) != n:
            raise ValueError(f"group_ids for '{varname}' has wrong length")

        group_ids_list.append(np.ascontiguousarray(gid, dtype=np.int32))
        cat_counts_list.append(ncat)
        targets_list.append(np.ascontiguousarray(props, dtype=np.float64))

    # Build initial weights
    if start_weights is not None:
        w = np.ascontiguousarray(start_weights, dtype=np.float64)
        w = w * len(w) / w.sum()
    else:
        w = np.ones(n, dtype=np.float64)

    params = {
        "min_weight":     min_weight,
        "max_weight":     max_weight,
        "inner_max_iter": max_iterations,
        "outer_max_iter": max_iterations,  # mirrors R bridge: outer = inner = max_iterations
        # legacy tol_abs: use absolute_tol when set, else fall back to 1e-6 for old C path
        "tol_abs":        absolute_tol if absolute_tol > 0.0 else 1e-6,
        "verbose":        verbose,
        "algorithm":      alg_int,
        "epsilon":        0.05,
        "bounds_mode":    _bounds_mode_int,
        # Convergence config (WU-G)
        "pct_tol":        pct_tol,
        "absolute_tol":   absolute_tol,
        "metric":         metric,    # replaces criterion
        "rule":           rule,      # new
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
    }
    if capacity_penalty is not None:
        params["capacity_penalty"] = capacity_penalty
    if alm_penalty is not None and alm_penalty > 0.0:
        params["alm_penalty"] = alm_penalty

    log_fn = print if verbose > 0 else None

    _, weights_out, result_dict = calibrate(
        n, K, w, group_ids_list, cat_counts_list, targets_list, params, log_fn
    )

    # Build convergence_used nested dict from raw integer fields (WU-G).
    # Mirrors R's harvest.R WU-E2 block.
    result_dict["convergence_used"] = {
        "metric":        _METRIC_NAMES[result_dict.get("convergence_metric", 0)],
        "rule":          _RULE_NAMES[result_dict.get("convergence_rule", 1)],
        "tol":           result_dict.get("convergence_tol", 0.001),
        "fired_at_iter": result_dict.get("convergence_iter", -1),
        "objective":         result_dict.get("convergence_solver_objective", float("inf")),
        "minimized_metric":  _METRIC_NAMES[result_dict.get("convergence_minimized_metric", 0)],
    }

    if result_dict["status"] == 1:
        warnings.warn(
            f"leafblower: calibration did not converge (max_error={result_dict['max_error']:.2e})",
            UserWarning, stacklevel=2
        )
    if result_dict.get("n_bounds_violated", 0) > 0:
        warnings.warn(
            f"cell-mode bounds: {result_dict['n_bounds_violated']} weights fell outside "
            f"[{min_weight:.3f}, {max_weight:.3f}]. Consider bounds_mode='unit'.",
            UserWarning, stacklevel=2)
    if result_dict.get("n_bounds_clamped", 0) > 0:
        warnings.warn(
            f"unit-mode bounds: {result_dict['n_bounds_clamped']} weights clamped "
            f"to [{min_weight:.3f}, {max_weight:.3f}] during per-cell water-filling.",
            UserWarning, stacklevel=2)
    if result_dict["status"] == 2:
        raise RuntimeError(
            "leafblower: infeasible problem — persistent empty cell with positive target "
            "(detected after 5 consecutive outer iterations)"
        )
    elif result_dict["status"] == 3:
        raise ValueError(f"leafblower: invalid arguments — {result_dict['message']}")
    if result_dict["status"] == 4:  # RK_ERR_BUDGET
        e_final   = result_dict.get("best_error", float("nan"))
        b_iter    = result_dict.get("best_iter", 0)
        iters     = result_dict.get("iterations", 0)
        mstr      = result_dict.get("convergence_used", {}).get("metric", "metric")
        tol_used  = absolute_tol if absolute_tol > 0.0 else pct_tol

        stall_ratio = (b_iter / iters) if iters > 0 else 1.0
        if stall_ratio < 0.5 and math.isfinite(e_final):
            warnings.warn(
                f"leafblower: fixed point at {mstr}={e_final:.2e} "
                f"(best at iter {b_iter} of {iters}, ratio={stall_ratio:.2f}). "
                f"More iterations will not improve calibration. "
                f"Try: accelerate=True, method='newton_kl', or method='ieppa+accel'.",
                UserWarning, stacklevel=2)
        else:
            e_prev    = result_dict.get("metric_prev_check", float("inf"))
            prev_iter = result_dict.get("prev_check_iter", -1)
            interval  = b_iter - prev_iter
            has_prev  = (math.isfinite(e_prev) and math.isfinite(e_final) and
                         e_prev > e_final > 0 and
                         interval > 0 and
                         math.isfinite(tol_used) and tol_used > 0)
            if has_prev:
                r_est = (e_final / e_prev) ** (1.0 / interval)
                if 0 < r_est < 1:
                    n_more  = math.ceil(math.log(tol_used / e_final) / math.log(r_est))
                    n_total = b_iter + n_more
                    warnings.warn(
                        f"leafblower: budget exhausted — {mstr}={e_final:.2e} "
                        f"at {iters} iters. Asymptotic rate r={r_est:.4f} "
                        f"(last {interval} iters): ~{n_total:.0f} total iterations needed.",
                        UserWarning, stacklevel=2)
                else:
                    warnings.warn(
                        f"leafblower: budget exhausted — {mstr}={e_final:.2e} "
                        f"at {iters} iters. Increase max_iterations.",
                        UserWarning, stacklevel=2)
            else:
                warnings.warn(
                    f"leafblower: budget exhausted — {mstr}={e_final:.2e} "
                    f"at {iters} iters. Increase max_iterations.",
                    UserWarning, stacklevel=2)

    # Solver returns sum(weights) = n (enforced in src/ieppa.cpp, src/raking.cpp,
    # src/lbfgsb_solver.cpp per user directive 2026-04-24). No wrapper-level
    # normalization — removing it preserves the bounds_mode="unit" strict-bounds
    # guarantee (ieppa's water-fill clamps are final; not re-pushed by post-scale).

    # NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping
    # here would break sum(weights * d) == target totals when per-cell mixing
    # parameters d are non-uniform: individual weights may legitimately exceed
    # per-cell bounds after expansion even when cell aggregates are in range.
    # The iEPPA/LBFGSB solvers enforce bounds on the cell aggregate X[c], which
    # is the invariant that preserves calibration. See
    # tests/testthat/test-ieppa-nonuniform-d.R.

    # weights_out is already a copy (contract from _bindings.cpp)
    if not attach_weights:
        weights_out_arr = np.asarray(weights_out, dtype=np.float64)
        # Expose result diagnostics even when attach_weights=False.
        # ndarray.attrs does not exist, so we must return a dict wrapper with
        # 'weights' key containing the array and 'result' key for diagnostics.
        # This preserves the contract: callers can access diagnostics via
        # result['result'] when attach_weights=False (same key as DataFrame.attrs).
        return {
            "weights": weights_out_arr,
            "result": result_dict,
        }

    if _PANDAS_AVAILABLE:
        out = data.copy()
        out[weight_column] = weights_out
        # Expose calibration diagnostics via DataFrame.attrs (PEP 526 / pandas 1.0+).
        # Nest SOR fields to match R's result$sor namespace.
        result_dict["sor"] = {
            "min_omega": result_dict.pop("sor_min_omega"),
            "n_damped":  result_dict.pop("sor_n_damped"),
        }
        out.attrs["result"] = result_dict
        out.attrs["iterations"] = result_dict["iterations"]
        return out
    return weights_out


def diagnose_weights(data, targets, weights):
    """
    Diagnose calibration quality (Python equivalent of R diagnose_weights()).

    Parameters
    ----------
    data : pd.DataFrame
    targets : dict of dicts, e.g. {"age": {"18-34": 0.3, "35+": 0.7}}
    weights : 1D array-like of calibrated weights, length len(data)

    Returns
    -------
    pd.DataFrame with columns:
        variable, level, prop_original, prop_weighted, target,
        error_original, error_weighted
    """
    if not _PANDAS_AVAILABLE:
        raise ImportError("pandas required for diagnose_weights; pip install pandas")

    weights = np.asarray(weights, dtype=float)
    if len(weights) != len(data):
        raise ValueError("weights length must equal len(data)")

    rows = []
    for varname, tgt_dict in targets.items():
        series = data[varname]
        # Use pd.isna for null detection — avoids coercing actual "nan" strings to NA
        not_na = ~pd.isna(series)
        col_str = series.astype(str)  # convert after NA check
        col_notna = not_na.values
        w_total = weights[col_notna].sum()
        n_total = col_notna.sum()
        for level, tgt_val in tgt_dict.items():
            mask = col_notna & (col_str.values == str(level))
            n_lvl = mask.sum()
            prop_orig = n_lvl / n_total if n_total > 0 else 0.0
            w_lvl = weights[mask].sum()
            prop_wtd = w_lvl / w_total if w_total > 0 else 0.0
            rows.append({
                "variable":       varname,
                "level":          level,
                "prop_original":  prop_orig,
                "prop_weighted":  prop_wtd,
                "target":         tgt_val,
                "error_original": prop_orig - tgt_val,
                "error_weighted": prop_wtd - tgt_val,
            })
    return pd.DataFrame(rows)
