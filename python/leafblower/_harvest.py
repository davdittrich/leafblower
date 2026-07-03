from __future__ import annotations
import math
import warnings
import numpy as np
import pandas as pd
from typing import Dict, Optional

from ._leafblower import calibrate

_TARGET_SUM_TOL = 5e-7  # matches R harvest() tolerance for target proportions summing to 1


def _validate_pos_scalar(name: str, val: float, *, null_clause: str = "",
                         warn_suffix: str = None, check_upper: bool = True) -> None:
    """Mirror R harvest.R positive-finite-scalar validation, per-param (R has 3
    distinct rules — capacity_penalty:378-393, alm_penalty:395-405,
    newton_tsvd_ratio:406-411).

    - Raises ValueError on non-finite / <=0 (all three); and on >1e15 when
      check_upper (capacity/alm only — newton_tsvd_ratio has no upper bound in R).
    - null_clause is R's per-param stop-message prefix: "NULL (auto) or " for
      capacity, "NULL (disabled) or " for alm, "" for newton_tsvd_ratio.
    - warn_suffix (when not None) mirrors R's per-param warning() for 0<val<1e-15;
      newton_tsvd_ratio has no such warning in R, so it passes warn_suffix=None.
    """
    if not math.isfinite(val) or val <= 0 or (check_upper and val > 1e15):
        raise ValueError(
            f"{name} must be {null_clause}a positive finite scalar; got: {val!r}"
        )
    if warn_suffix is not None and val < 1e-15:
        warnings.warn(
            f"{name}={val} is below recommended range; {warn_suffix}",
            UserWarning, stacklevel=2)

_METRIC_MAP = {
    "max_err": 0, "mean_err": 1, "kl": 2, "chi2": 3,
    "grake_norm": 4, "l1_weight": 5, "marginal_kl": 6,
}
_RULE_MAP = {"threshold": 0, "improvement": 1, "plateau": 2}
_STOP_WHEN_MAP = {"any": 0, "all": 1}
_KNOWN_CONVERGENCE_KEYS = frozenset({
    "metric", "rule", "tol", "pct", "absolute", "improvement", "stop_when",
})
_KNOWN_SOR_KEYS = frozenset({
    "auto", "omega_min", "omega_max", "omega", "omega_init", "burnin",
    "omega_mode_id",
})
_OMEGA_MODE_MAP = {"heuristic": 0, "fixed": 1, "spectral": 2}
# Reverse-lookup arrays mirror CalibMetric/CalibRule enum order in leafblower.h
_METRIC_NAMES = ["max_err", "mean_err", "kl", "chi2", "grake_norm", "l1_weight", "marginal_kl"]
_RULE_NAMES   = ["threshold", "improvement", "plateau"]


def _compute_sparseness_diag(df, targets, cat_threshold=0.01, obs_threshold=30,
                             na_injected=None):
    """Flag categories where T_kj < cat_threshold or n_kj < obs_threshold.

    na_injected: set of variable names that received an 'NA' bin from
    add_na_proportion=True. For those variables, the 'NA' level is counted
    via pd.isna() because value_counts(dropna=True) drops NaN/None entirely.
    """
    na_injected = na_injected or set()
    sparse_cats = {}
    for v, tgt in targets.items():
        if v not in df.columns:
            continue
        counts = df[v].astype(str).value_counts(dropna=True)
        for level, T_kj in tgt.items():
            lv_str = str(level)
            if lv_str == "NA" and v in na_injected:
                # NA-bin injected by add_na_proportion: value_counts(dropna=True)
                # drops NaN/None so they never appear under any string key; count
                # directly. Conflate true-NA + literal-string-"NA" into the "NA"
                # level to mirror the solver encoding (_harvest.py group_ids
                # else-branch: col.astype(str).where(~isna, other="NA")). The added
                # literal term is 0 when no row holds the string "NA".
                n_kj = int(pd.isna(df[v]).sum()) + int(
                    ((~pd.isna(df[v])) & (df[v].astype(str) == "NA")).sum()
                )
            else:
                n_kj = int(counts.get(lv_str, 0))
            if T_kj < cat_threshold or n_kj < obs_threshold:
                sparse_cats.setdefault(v, []).append(
                    {"level": level, "T_kj": T_kj, "n_kj": n_kj}
                )
    return sparse_cats


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
    """sor=None disables SOR (matches R default sor=NULL). sor=False/sor={} = explicit disable."""
    if sor is None:
        return {}  # R default: sor=NULL (disabled)
    if sor is False:
        return {}  # explicit disable
    return sor  # passthrough dict


def _parse_sor(sor):
    """Mirror R parse_sor(): returns (enabled, auto, omega_init, omega_min, omega_max, omega_fixed, burnin, omega_mode_id)."""
    if not sor:  # empty dict = disabled (from _resolve_sor(None) or _resolve_sor(False))
        return 0, 0, 1.0, 0.3, 1.5, -1.0, 20, None   # None = let C init default govern
    # bad keys in user-input order; valid keys in R's declared order (harvest.R:943-949)
    unknown = [k for k in sor if k not in _KNOWN_SOR_KEYS]
    if unknown:
        raise ValueError(
            f"Unknown sor key(s): {', '.join(unknown)}. "
            f"Valid keys: auto, omega_min, omega_max, omega, omega_init, burnin, omega_mode_id")
    enabled = 1
    auto = 1 if sor.get("auto", True) else 0
    omega_init = float(sor.get("omega_init", 1.0))
    omega_min = float(sor.get("omega_min", 0.3))
    omega_max = float(sor.get("omega_max", 1.5))
    omega_fixed = float(sor.get("omega", -1.0))
    burnin = int(sor.get("burnin", 20))
    omega_mode_id = sor.get("omega_mode_id", None)  # None = C rk_params_init default governs
    if omega_mode_id is not None:
        if isinstance(omega_mode_id, str):
            if omega_mode_id not in _OMEGA_MODE_MAP:
                raise ValueError(
                    f"Unknown omega_mode_id string '{omega_mode_id}'. "
                    f"Use 'heuristic', 'fixed', or 'spectral'.")
            omega_mode_id = _OMEGA_MODE_MAP[omega_mode_id]
        else:
            omega_mode_id = int(omega_mode_id)
    return enabled, auto, omega_init, omega_min, omega_max, omega_fixed, burnin, omega_mode_id


def _parse_target(target, target_map=None):
    """Mirror R parse_target(): convert DataFrame target to dict. No normalization."""
    if not isinstance(target, pd.DataFrame):
        return target  # already dict
    if target_map is not None:
        vcol = target_map["variable"]
        lcol = target_map["level"]
        pcol = target_map["proportion"]
    elif all(c in target.columns for c in ("variable", "level", "proportion")):
        vcol, lcol, pcol = "variable", "level", "proportion"
    elif target.shape[1] == 3:
        raise ValueError(
            "target DataFrame has 3 columns but no 'variable'/'level'/'proportion' names. "
            "Add column names or pass target_map=dict(variable=..., level=..., proportion=...).")
    else:
        raise ValueError(
            "Cannot determine variable/level/proportion columns in target DataFrame.")
    result = {}
    for v in target[vcol].unique():
        sub = target[target[vcol] == v]
        # CR-E13 (5ye4.13): a duplicated (variable, level) row is ambiguous input;
        # dict(zip(...)) would silently keep only the last. Reject in both R + Python.
        # Use pandas duplicated() (NaN==NaN, matching R's anyDuplicated) — a plain
        # set() would treat distinct np.nan objects as unique and miss NaN-level dups.
        dup_mask = sub[lcol].duplicated(keep=False)
        if dup_mask.any():
            dups = sorted(sub[lcol][dup_mask].unique(), key=str)
            raise ValueError(
                f"target for '{v}' has duplicate level(s): {dups}")
        d = dict(zip(sub[lcol], sub[pcol]))
        total = sum(d.values())
        if abs(total - 1.0) > _TARGET_SUM_TOL:
            raise ValueError(
                f"target proportions for '{v}' sum to {total:.6f}, not 1.0")
        result[v] = d
    return result


def harvest(
    data,
    targets: Dict[str, Dict[str, float]],
    min_weight: float = 0.0,
    max_weight: float = 5.0,
    method: str = "oris",
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
    target_map=None,
    ridge_lambda: float = 0.0,
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
    method : str, one of "oris" (default), "raking", "sinkhorn",
        "chebyshev", "greg", "oris_soft", "greenkhorn", "logit", "newton_kl"
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
        Note: method="chebyshev" ignores "rule" (and "metric"/"stop_when"); its
        interior-point solver stops on its own complementarity-gap criterion,
        falling back to the "pct"/"absolute" tolerance on the max marginal error.
    sor : dict for SOR adaptive under-relaxation (ORIS and raking). None disables SOR. Keys:
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
    # Removed autumn legacy params — raise TypeError with migration hint
    _REMOVED_PARAMS = {
        "select_params":  "use method= and scheduler= instead",
        "select_function": "not supported in leafblower; use method=",
        "error_function":  "not supported in leafblower",
        "adaptive_order":  "use scheduler='greedy' instead",
        "enforce_mean":    "leafblower always enforces sum(w)=n; param removed",
    }
    for _p, _hint in _REMOVED_PARAMS.items():
        if _p in _kwargs:
            raise TypeError(f"harvest() got removed autumn param '{_p}': {_hint}")

    # CR-E4 (5ye4.4): reject any remaining unrecognized kwargs. The catch-all
    # **_kwargs previously swallowed typos (max_weigth=3.0, bounds_mod="unit")
    # silently, producing wrong results with no signal the argument was ignored.
    # Known forward-compat names are the _REMOVED_PARAMS above (raised already);
    # anything left is a caller mistake.
    if _kwargs:
        raise TypeError(
            f"harvest() got unexpected keyword argument(s): {sorted(_kwargs)}"
        )

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
    sor_enabled, sor_auto, sor_omega_init, sor_omega_min, sor_omega_max, sor_omega_fixed, sor_burnin, sor_omega_mode_id = _parse_sor(
        _resolve_sor(sor)
    )

    # Convert dict data to DataFrame; validate input type.
    if isinstance(data, dict):
        data = pd.DataFrame(data)
    elif not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pd.DataFrame or dict")

    n = len(data)
    if n == 0:  # mirror R harvest.R:420-421
        raise ValueError("leafblower: 'data' must be a non-empty data.frame")

    # Method mapping
    method_lc = method.lower()

    # Mirror R harvest.R:308-331: per-method natural convergence objective.
    # Only fires when user gave no explicit metric/improvement/pct/absolute.
    _conv = convergence or {}
    if not any(k in _conv for k in ("metric", "improvement", "pct", "absolute")):
        _method_metric = {
            "oris": "marginal_kl", "oris_soft": "marginal_kl",
            "raking": "kl", "greenkhorn": "kl",
            "sinkhorn": "kl", "newton_kl": "kl",
            "greg": "chi2",
            # chebyshev: omitted — max_err is its natural objective (L-inf)
            # logit: omitted — no natural KL objective, keep max_err
        }.get(method_lc)
        if _method_metric is not None:
            metric = _METRIC_MAP[_method_metric]

    alg_map = {
        "oris": 1, "raking": 3,
        "sinkhorn": 4, "chebyshev": 5, "greg": 6,
        "oris_soft": 8, "greenkhorn": 9, "logit": 10, "newton_kl": 11,
    }  # "auto" (0) removed from Python user API; "grake" (7) removed; "lbfgsb" (2) removed — enum gap
    if method_lc not in alg_map:
        raise ValueError(f"method must be one of {sorted(alg_map)}; 'auto' is R-only and not supported in Python")
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

    # CR-C15 (kxna.15): reject max_iterations < 1, mirroring R harvest.R. R's wrapper
    # already guards this; Python bypassed it, so max_iterations=0 collapsed the solver
    # loop to zero iterations and (for logit) returned all-zero weights with STALL — a
    # silently wrong result rather than an error.
    if not isinstance(max_iterations, (int, np.integer)) or max_iterations < 1:
        raise ValueError(
            f"max_iterations must be a positive integer; got: {max_iterations!r}"
        )

    # Convert DataFrame target to dict (mirrors R parse_target())
    targets = _parse_target(targets, target_map)

    # 81bx: auto_collapse — merge rare categories into __other__.
    # CR-E5 (5ye4.5): runs AFTER dict→DataFrame coercion and _parse_target so
    # data[v].astype(str) and targets[v].items() operate on real DataFrame
    # columns / parsed dict-of-floats (was crashing on dict-data / DataFrame-targets).
    # CR-E7 (5ye4.7): bool(auto_collapse) — truthy 1 enables, 0/False/None disable
    # (was `is True`, which silently ignored auto_collapse=1).
    if bool(auto_collapse) or (collapse_vars is not None and auto_collapse is not False):
        vars_to_collapse = collapse_vars if collapse_vars is not None else list(targets.keys())
        data = data.copy()  # avoid mutating caller's DataFrame
        # Deep-copy targets too: targets[v].pop() below permanently corrupts
        # caller's dict (parity with the data.copy() above).
        targets = {k: (dict(v) if isinstance(v, dict) else v) for k, v in targets.items()}
        for v in vars_to_collapse:
            if v not in targets:
                continue
            # CR-E6 (5ye4.6): hoist the O(n) column scan out of the per-level
            # loop — one value_counts() per variable, reused across levels,
            # instead of a full `astype(str) == lv` rescan per level (was
            # O(n × levels)). Matches the _compute_sparseness_diag pattern.
            # dropna=True is irrelevant here: astype(str) already turns NaN into
            # the "nan" string, and "NA" levels are excluded from `rare` below,
            # so no separate pd.isna() count is needed for rare determination.
            counts = data[v].astype(str).value_counts(dropna=True)
            rare = [
                lv for lv, t_kj in targets[v].items()
                if (t_kj < 0.01 or int(counts.get(str(lv), 0)) < 30)
                and str(lv) != "NA"
            ]
            if not rare:
                continue
            other_mass = sum(targets[v].pop(lv) for lv in rare)
            targets[v]["__other__"] = targets[v].get("__other__", 0.0) + other_mass
            rare_set = {str(lv) for lv in rare}
            col = data[v].astype(str)
            data[v] = col.where(~col.isin(rare_set), other="__other__")

    # yaye: add_na_proportion — renormalize targets and add explicit NA bin
    _na_vars: set = set()
    # CR-E7 (5ye4.7): bool() — 0/False/None disable NA-bin injection (was
    # `is not False`, which treated 0 and None as enabled, the opposite of intent).
    if bool(add_na_proportion):
        for v in list(targets.keys()):
            if v not in data.columns:
                continue
            na_frac = float(pd.isna(data[v]).mean())
            if na_frac == 0.0:
                continue
            if na_frac == 1.0:
                raise ValueError(
                    f"add_na_proportion: all observations are NA for margin '{v}'"
                )
            targets[v] = {k: val * (1.0 - na_frac) for k, val in targets[v].items()}
            targets[v]["NA"] = na_frac
            _na_vars.add(v)

    # c8w1: sparseness diagnostic (pre-solve)
    _sparse_diag = _compute_sparseness_diag(data, targets, cat_threshold=0.01, obs_threshold=30,
                                             na_injected=_na_vars)
    if _sparse_diag:
        n_flagged = sum(len(v) for v in _sparse_diag.values())
        warnings.warn(
            f"leafblower: {n_flagged} sparse categories detected "
            f"(T_kj < 1% or n_kj < 30); see result['diagnostics']['sparseness']",
            UserWarning, stacklevel=2,
        )

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

        # Vectorized encoding via pandas Categorical: O(n) in C, ~6x faster than
        # a Python for-loop. col.astype(str) handles mixed-type columns; pd.isna
        # entries map to codes=-1 automatically when observed=False.
        # yaye: for NA-bin margins, fill NAs with "NA" so they map to the NA bin.
        if varname in _na_vars:
            col_enc = col.astype(str).where(~pd.isna(col), other="NA")
        else:
            col_enc = col.astype(str).where(~pd.isna(col), other=np.nan)
        cat = pd.Categorical(col_enc, categories=[str(lv) for lv in levels])
        gid = cat.codes.astype(np.int32)  # -1 for NA/unknown levels

        if len(gid) != n:
            raise ValueError(f"group_ids for '{varname}' has wrong length")

        group_ids_list.append(np.ascontiguousarray(gid, dtype=np.int32))
        cat_counts_list.append(ncat)
        targets_list.append(np.ascontiguousarray(props, dtype=np.float64))

    # Build initial weights
    if start_weights is not None:
        w = np.ascontiguousarray(start_weights, dtype=np.float64)
        # CR-E10 (5ye4.10): validate BEFORE the rescale. `w * len(w) / w.sum()`
        # silently yields all-NaN when w.sum()==0, and a mis-length / negative /
        # non-finite start_weights otherwise surfaces as a confusing error far
        # from the bad-input site. Fail fast, naming the parameter.
        if w.shape[0] != n:
            raise ValueError(
                f"start_weights has length {w.shape[0]}, expected {n} "
                "(one per observation)"
            )
        if not np.isfinite(w).all():
            raise ValueError("start_weights contains non-finite values (NaN or inf)")
        if (w < 0).any():
            raise ValueError("start_weights contains negative values")
        w_sum = w.sum()
        if w_sum <= 0.0:
            raise ValueError("start_weights sums to zero or less; cannot rescale to n")
        w = w * len(w) / w_sum
    else:
        w = np.ones(n, dtype=np.float64)

    # Mirror R harvest.R:414-417: accelerate only supported for these methods.
    _ACCELERATE_METHODS = {"raking", "greenkhorn", "oris", "oris_soft"}
    if accelerate and method_lc not in _ACCELERATE_METHODS:
        warnings.warn(
            f"accelerate=True is only supported for method='raking', 'greenkhorn', "
            f"'oris', or 'oris_soft'; ignoring for method='{method}'",
            UserWarning, stacklevel=2)
        accelerate = False

    # Positive-finite-scalar validation, per-param to match R exactly.
    # capacity_penalty / alm_penalty are validated only when not None (their
    # disable path); an explicit -1.0 is rejected (<=0) just like R. Each mirrors
    # R's distinct message/rules: capacity (harvest.R:378-393, "NULL (auto)"),
    # alm (395-405, "NULL (disabled)", "objective penalty"), newton_tsvd_ratio
    # (406-411, unconditional, NOT nullable, no upper bound, no small-value warn).
    if capacity_penalty is not None:
        _validate_pos_scalar("capacity_penalty", capacity_penalty,
                             null_clause="NULL (auto) or ",
                             warn_suffix="constraint enforcement may be ineffective")
    if alm_penalty is not None:
        _validate_pos_scalar("alm_penalty", alm_penalty,
                             null_clause="NULL (disabled) or ",
                             warn_suffix="objective penalty may be ineffective")
    _validate_pos_scalar("newton_tsvd_ratio", newton_tsvd_ratio,
                         null_clause="", warn_suffix=None, check_upper=False)

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
        "sor_omega_max":  sor_omega_max,
        "sor_omega_fixed": sor_omega_fixed,
        "sor_burnin":     sor_burnin,
        "sor_omega_mode_id": sor_omega_mode_id,  # None = omit (C default governs)
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
    if capacity_penalty is not None:
        params["capacity_penalty"] = capacity_penalty
    if alm_penalty is not None:
        # Pass through unconditionally — R sends -1.0 as 'disabled' sentinel;
        # C-side handles sentinel semantics. Previously the >0.0 gate silently
        # diverged from R parity for negative inputs.
        params["alm_penalty"] = float(alm_penalty)

    log_fn = print if verbose > 0 else None

    _, weights_out, result_dict = calibrate(
        n, K, w, group_ids_list, cat_counts_list, targets_list, params, log_fn
    )

    # Build convergence_used nested dict from raw integer fields (WU-G).
    # Mirrors R's harvest.R WU-E2 block.
    def _safe_metric_name(idx):
        if 0 <= idx < len(_METRIC_NAMES):
            return _METRIC_NAMES[idx]
        raise RuntimeError(
            f"C bridge returned out-of-range metric id {idx} "
            f"(valid: 0..{len(_METRIC_NAMES)-1}). Check leafblower build version."
        )
    result_dict["convergence_used"] = {
        "metric":        _safe_metric_name(result_dict.get("convergence_metric", 0)),
        "rule":          _RULE_NAMES[result_dict.get("convergence_rule", 1)],
        "tol":           result_dict.get("convergence_tol", 0.001),
        "fired_at_iter": result_dict.get("convergence_iter", -1),
        "objective":         result_dict.get("convergence_solver_objective", float("inf")),
        "minimized_metric":  _safe_metric_name(result_dict.get("convergence_minimized_metric", 0)),
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
        # eb79.23: surface the solver's own message (naming the infeasible margin, e.g.
        # logit's structural pre-check), matching R harvest.R:645 and the status==3 branch
        # below. Fall back to a generic string only when the solver set no message (sinkhorn/
        # oris currently leave it empty on INFEAS — behavior unchanged for them).
        msg = result_dict.get("message") or ""
        raise RuntimeError(
            f"leafblower: {msg}" if msg else "leafblower: infeasible problem"
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
                f"Try: accelerate=True, method='newton_kl', or method='oris' with accelerate=True.",
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
    if result_dict["status"] == 5:  # RK_ERR_STALL
        # CR-E3 (5ye4.3): surface STALL (was silent — only statuses 1-4 branched).
        # Mirrors R harvest.R:729-734: constrained-optimum plateau; weights valid.
        # accelerate=True → SRAA-m weight-change plateau variant; else loss-function
        # plateau variant. (R checks the same `accelerate` flag.)
        if accelerate:
            warnings.warn(
                "leafblower: SRAA-m weight-change plateau — at constrained optimum; "
                "weights are valid; no further improvement is achievable",
                UserWarning, stacklevel=2)
        else:
            warnings.warn(
                "leafblower: loss function plateau — at constrained optimum given bounds; "
                "weights are valid; no further improvement is achievable",
                UserWarning, stacklevel=2)

    # Solver returns sum(weights) = n (enforced at exit by
    # lbw::finalize_weights / lbw::finalize_weights_buf in src/calib_dispatch.hpp).
    # No wrapper-level
    # normalization — removing it preserves the bounds_mode="unit" strict-bounds
    # guarantee (oris's water-fill clamps are final; not re-pushed by post-scale).

    # NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping
    # here would break sum(weights * d) == target totals when per-cell mixing
    # parameters d are non-uniform: individual weights may legitimately exceed
    # per-cell bounds after expansion even when cell aggregates are in range.
    # The ORIS solver enforces bounds on the cell aggregate X[c], which
    # is the invariant that preserves calibration. See
    # tests/testthat/test-oris-nonuniform-d.R.

    # c8w1: attach sparseness diagnostics to result_dict
    result_dict["diagnostics"] = {
        "sparseness": {
            "sparse_categories":  _sparse_diag,
            "pct_bounds_clamped": (result_dict.get("n_bounds_clamped", 0) / n
                                   if n > 0 else 0.0),
            "thresholds": {"target": 0.01, "obs": 30},
        }
    }

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

    out = data.copy()
    out[weight_column] = weights_out
    # Expose calibration diagnostics via DataFrame.attrs (PEP 526 / pandas 1.0+).
    # Nest SOR fields to match R's result$sor namespace.
    # Non-SOR solvers (newton_kl, chebyshev, greg, etc.) don't populate
    # SOR fields; .pop with default avoids KeyError. R side returns NULL
    # for these; Python returns None inside the nested dict.
    result_dict["sor"] = {
        "min_omega":    result_dict.pop("sor_min_omega", None),
        "n_damped":     result_dict.pop("sor_n_damped", None),
        "omega_mean":   result_dict.pop("sor_omega_mean", None),
        "any_latched":  result_dict.pop("sor_any_latched", None),
        "n_pinned_fb":  result_dict.pop("sor_n_pinned_fb", None),
        "n_warmup_fb":  result_dict.pop("sor_n_warmup_fb", None),
        "n_conv_fb":    result_dict.pop("sor_n_conv_fb", None),
        "n_resid_grew": result_dict.pop("sor_n_resid_grew", None),
        "n_monotone_cd": result_dict.pop("sor_n_monotone_cd", None),
    }
    out.attrs["result"] = result_dict
    out.attrs["iterations"] = result_dict["iterations"]
    return out


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

    Notes
    -----
    An "NA" entry in a target is treated as the missing-data bin: a row counts
    toward it if it is NA or its value is the literal string "NA" (conflation),
    matching harvest(add_na_proportion=True). A hand-built target naming a real
    "NA" category therefore also counts genuinely-missing rows in that bin.
    """
    weights = np.asarray(weights, dtype=float)
    if len(weights) != len(data):
        raise ValueError("weights length must equal len(data)")

    rows = []
    for varname, tgt_dict in targets.items():
        series = data[varname]
        # Use pd.isna for null detection — avoids coercing actual "nan" strings to NA
        na_mask = pd.isna(series).values
        col_str = series.astype(str)  # convert after NA check
        has_na_bin = "NA" in {str(lv) for lv in tgt_dict}
        if has_na_bin:
            # add_na_proportion case: all-obs denominators. NA observations are
            # counted via the explicit "NA" bin so shares sum to 1 when every obs
            # falls in a named level or the NA bin (an out-of-vocabulary value
            # lands in no bin -> sum < 1).
            w_total = weights.sum()
            n_total = len(series)
        else:
            # No "NA" bin (common case): exclude NA observations from BOTH the
            # level masks and the denominators. harvest drops NA/gid<0 obs from
            # the marginal constraints (raking.cpp: `if (g>=0)`), so named-level
            # shares must be measured over non-NA obs only — otherwise they sum
            # to <1 and produce a spurious error_weighted ~ -na_frac*target on
            # every level for well-calibrated data.
            not_na = ~na_mask
            w_total = weights[not_na].sum()
            n_total = int(not_na.sum())
        for level, tgt_val in tgt_dict.items():
            # The injected add_na_proportion bin (literally named "NA") CONFLATES
            # true-missing rows with rows whose value is the literal string "NA"
            # (na_mask OR col_str=="NA"), matching the solver encoding
            # (R harvest.R:130-131,475; _harvest.py:404-413). A literal-"NA" row
            # falls into the NA bin, not its own level. Supersedes 4ihf.4
            # (mask-only), which under-reported the NA bin and broke Σshares==1.
            if str(level) == "NA":
                mask = na_mask | (col_str.values == "NA")
            else:
                mask = (~na_mask) & (col_str.values == str(level))
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
