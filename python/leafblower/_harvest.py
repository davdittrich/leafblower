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


def harvest(
    data,
    targets: Dict[str, Dict[str, float]],
    min_weight: float = 0.0,
    max_weight: float = 5.0,
    method: str = "auto",
    verbose: int = 0,
    max_iterations: int = 500,
    start_weights: Optional[np.ndarray] = None,
    attach_weights: bool = True,
    weight_column: str = "weights",
    convergence: Optional[Dict] = None,
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
    method : "auto" | "ieppa" | "lbfgsb"
    verbose : int, 0=silent, 1=progress, 2=debug
    max_iterations : int, inner BCD max sweeps per outer iter (default 500)
    start_weights : optional 1D float64 array of initial weights
    attach_weights : if True, return DataFrame with weights column appended
    weight_column : name of the weights column (default "weights")
    convergence : dict; key "absolute" → tol_abs; "pct" → DeprecationWarning;
                  unknown keys → ValueError
    Returns
    -------
    pd.DataFrame (if attach_weights=True) or np.ndarray
    """
    # Handle convergence dict
    tol_abs = 1e-6
    if convergence is not None:
        for k in convergence:
            if k == "absolute":
                tol_abs = float(convergence[k])
            elif k == "pct":
                warnings.warn(
                    "convergence['pct'] is deprecated; use 'absolute' instead.",
                    DeprecationWarning, stacklevel=2,
                )
            else:
                raise ValueError(f"unknown convergence key '{k}'")

    # Convert dict data to DataFrame
    if isinstance(data, dict):
        if not _PANDAS_AVAILABLE:
            raise ImportError("pandas required to use dict input; install with pip install pandas")
        data = pd.DataFrame(data)

    if _PANDAS_AVAILABLE and not isinstance(data, pd.DataFrame):
        raise TypeError("data must be a pd.DataFrame or dict")

    n = len(data)

    # Method mapping
    method_lc = method.lower()
    if method_lc in ("rake", "nrake"):
        warnings.warn(f"method='{method_lc}' (IPF) not implemented; using L-BFGS-B", UserWarning, stacklevel=2)
        method_lc = "lbfgsb"
    elif method_lc == "nr":
        warnings.warn("method='nr' not implemented; using L-BFGS-B", UserWarning, stacklevel=2)
        method_lc = "lbfgsb"

    alg_map = {"auto": 0, "ieppa": 1, "lbfgsb": 2}
    if method_lc not in alg_map:
        raise ValueError(f"method must be one of {list(alg_map)}")
    alg_int = alg_map[method_lc]

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
        "tol_abs":        tol_abs,
        "verbose":        verbose,
        "algorithm":      alg_int,
        "epsilon":        0.05,
        "lbfgs_m":        10,
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
    elif result_dict["status"] == 2:
        raise RuntimeError("leafblower: infeasible problem — empty cell with positive target")
    elif result_dict["status"] == 3:
        raise ValueError(f"leafblower: invalid arguments — {result_dict['message']}")

    # weights_out is already a copy (contract from _bindings.cpp)
    if not attach_weights:
        return weights_out

    if _PANDAS_AVAILABLE:
        out = data.copy()
        out[weight_column] = weights_out
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
