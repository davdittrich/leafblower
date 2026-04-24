from __future__ import annotations
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

_KNOWN_CONVERGENCE_KEYS = frozenset({"pct", "absolute", "criterion", "stop_when"})
_CRITERION_MAP = {"pct": 0, "max_err": 1, "mean_err": 2, "kl": 3, "chi2": 4}
_STOP_WHEN_MAP = {"any": 0, "all": 1}


def _parse_convergence(conv):
    """Mirror R parse_convergence(): derive pct_tol, absolute_tol, criterion, stop_when."""
    if conv is None:
        conv = {}
    unknown = set(conv) - _KNOWN_CONVERGENCE_KEYS
    if unknown:
        raise ValueError(f"unknown convergence key '{next(iter(unknown))}'")
    explicit_pct = "pct" in conv
    explicit_abs = "absolute" in conv
    pct_tol = (
        float(conv["pct"]) if explicit_pct
        else (0.001 if not explicit_abs else 0.0)
    )
    absolute_tol = float(conv.get("absolute", 0.0))
    default_criterion = "pct" if (explicit_pct or not explicit_abs) else "max_err"
    criterion_str = conv.get("criterion", default_criterion)
    if criterion_str not in _CRITERION_MAP:
        raise ValueError(f"convergence criterion must be one of {list(_CRITERION_MAP)}")
    criterion = _CRITERION_MAP[criterion_str]
    stop_when_str = conv.get("stop_when", "any")
    if stop_when_str not in _STOP_WHEN_MAP:
        raise ValueError(f"stop_when must be 'any' or 'all'")
    stop_when = _STOP_WHEN_MAP[stop_when_str]
    return pct_tol, absolute_tol, criterion, stop_when


def _parse_sor(sor):
    """Mirror R parse_sor(): returns (enabled, auto, omega_init, omega_min, omega_fixed, burnin)."""
    if sor is None:
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
    **kwargs,
):
    """
    Calibrate survey weights. Drop-in for R leafblower::harvest().

    Parameters
    ----------
    data : pd.DataFrame or dict of lists
    targets : dict of dicts, e.g. {"age": {"18-34": 0.3, "35+": 0.7}}
    min_weight : float, lower bound on weights (default 0 = no bound)
    max_weight : float, upper bound on weights (default 5)
    method : "ieppa" | "lbfgsb" (default "ieppa")
    verbose : int, 0=silent, 1=progress, 2=debug
    max_iterations : int, inner BCD max sweeps per outer iter (default 500)
    start_weights : optional 1D float64 array of initial weights
    attach_weights : if True, return DataFrame with weights column appended
    weight_column : name of the weights column (default "weights")
    convergence : dict controlling the stopping criterion. Keys:
        "pct" (float, default 0.001) — max proportional weight-change threshold.
        "absolute" (float) — absolute threshold for the active criterion.
        "criterion" (str) — one of "pct" (default), "max_err", "mean_err", "kl", "chi2".
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
    # Parse convergence and SOR
    pct_tol, absolute_tol, criterion, stop_when = _parse_convergence(convergence)
    sor_enabled, sor_auto, sor_omega_init, sor_omega_min, sor_omega_fixed, sor_burnin = _parse_sor(sor)

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
    if method_lc in ("rake", "nrake"):
        warnings.warn(f"method='{method_lc}' (IPF) not implemented; using L-BFGS-B", UserWarning, stacklevel=2)
        method_lc = "lbfgsb"
    elif method_lc == "nr":
        warnings.warn("method='nr' not implemented; using L-BFGS-B", UserWarning, stacklevel=2)
        method_lc = "lbfgsb"

    alg_map = {"ieppa": 1, "lbfgsb": 2, "raking": 3}  # "auto" (0) removed from user API
    if method_lc not in alg_map:
        raise ValueError(f"method must be one of {list(alg_map)}")
    alg_int = alg_map[method_lc]

    if bounds_mode not in ("cell", "unit"):
        raise ValueError(f"bounds_mode must be 'cell' or 'unit', got {bounds_mode!r}")
    _bounds_mode_int = {"cell": 0, "unit": 1}[bounds_mode]

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
        if abs(sum(props) - 1.0) > 1e-8:
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
        "lbfgs_m":        10,
        "bounds_mode":    _bounds_mode_int,
        # Convergence config (WU-F)
        "pct_tol":        pct_tol,
        "absolute_tol":   absolute_tol,
        "criterion":      criterion,
        "stop_when":      stop_when,
        # SOR config (WU-F)
        "sor_enabled":    sor_enabled,
        "sor_auto":       sor_auto,
        "sor_omega_init": sor_omega_init,
        "sor_omega_min":  sor_omega_min,
        "sor_omega_fixed": sor_omega_fixed,
        "sor_burnin":     sor_burnin,
    }

    log_fn = print if verbose > 0 else None

    _, weights_out, result_dict = calibrate(
        n, K, w, group_ids_list, cat_counts_list, targets_list, params, log_fn
    )

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
        weights_out_arr = np.array(weights_out)
        # Attach result metadata via a wrapper object attribute trick is not possible
        # for bare ndarray; return as-is (result_dict accessible via the tuple form
        # only when calling the low-level calibrate() directly).
        return weights_out

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
