"""Henry & Valliant (2015) calibration design effect — Python wrapper.

All computation lives in the C++17 core (src/design_effect.cpp).
This module marshals pandas/numpy data into the C ABI shape and returns a float.

Henry, K.A. & Valliant, R. (2015). A design effect measure for calibration
weighting in single-stage samples. Survey Methodology, 41(2), 315-331.
Statistics Canada Catalogue No. 12-001-X. Equation 3.5.
"""

from __future__ import annotations

import warnings
from typing import Mapping, Optional, TYPE_CHECKING

import numpy as np

from leafblower._leafblower import _design_effect as _c_design_effect

if TYPE_CHECKING:  # annotation-only; runtime pandas import is lazy (see design_effect body)
    import pandas as pd


def design_effect(
    weights: np.ndarray,
    outcome: Optional[np.ndarray] = None,
    data: Optional[pd.DataFrame] = None,
    target: Optional[Mapping[str, Mapping[str, float]]] = None,
) -> float:
    """Kish (1965) or Henry & Valliant (2015) Eq 3.5 design effect.

    Parameters mirror the R wrapper. See:
    - Kish, L. (1965). Survey Sampling. Wiley.
    - Henry, K.A. & Valliant, R. (2015). A design effect measure for
      calibration weighting in single-stage samples. Survey Methodology,
      41(2), 315-331. Statistics Canada Catalogue No. 12-001-X. Eq 3.5.

    All numerical computation is in C++17; R + Python wrappers produce
    byte-identical output by construction.
    """
    weights = np.ascontiguousarray(weights, dtype=np.float64)
    # CR-F7b (5ye4.16): enforce the documented weight contract (no NA, all finite,
    # sum > 0) on BOTH the 1-arg and 4-arg paths, mirroring R design_effect.R:41-46.
    # Previously neither path validated weights, silently passing NaN/Inf/zero-sum to C.
    if np.isnan(weights).any():
        raise ValueError("design_effect: weights must not contain NA values")
    if not np.isfinite(weights).all():
        raise ValueError("design_effect: weights must all be finite")
    if weights.sum() <= 0:
        raise ValueError("design_effect: sum(weights) must be positive")
    if outcome is None:
        return float(_c_design_effect(weights, None, None, None, 0)["deff_K"])
    if data is None or target is None:
        raise ValueError("design_effect: 4-argument form requires both 'data' and 'target'")
    n = weights.shape[0]
    outcome = np.ascontiguousarray(outcome, dtype=np.float64)
    if outcome.shape[0] != n:
        raise ValueError(
            f"design_effect: outcome length ({outcome.shape[0]}) != len(weights) ({n})"
        )
    if data.shape[0] != n:
        raise ValueError(
            f"design_effect: nrow(data) ({data.shape[0]}) != len(weights) ({n})"
        )
    K = len(target)
    if K == 0:
        return float(_c_design_effect(weights, outcome, None, None, 0)["deff_K"])
    missing = [v for v in target if v not in data.columns]
    if missing:
        raise ValueError(f"design_effect: data missing target column(s): {missing}")
    try:
        import pandas as pd
    except ImportError as exc:
        raise ImportError(
            "design_effect with 'data'/'target' arguments requires pandas. "
            "Install it with: pip install pandas"
        ) from exc
    data_codes = np.empty(n * K, dtype=np.int32)
    cat_counts = np.empty(K, dtype=np.int32)
    for k, var in enumerate(target):
        levs = list(target[var].keys())
        col = data[var]
        if col.isna().any():
            raise ValueError(f"design_effect: data[{var!r}] contains NA")
        bad = set(col.unique()) - set(levs)
        if bad:
            raise ValueError(
                f"design_effect: data[{var!r}] has level(s) {sorted(bad)} not in target"
            )
        cat = pd.Categorical(col, categories=levs)
        codes = np.array(cat.codes, dtype=np.int32)
        for i in range(n):
            data_codes[i * K + k] = codes[i]
        cat_counts[k] = len(levs)
    res = _c_design_effect(weights, outcome, data_codes, cat_counts, K)
    if res["rank_def"] == 1:
        warnings.warn(
            "design_effect: calibration margins rank-deficient; deff_H = deff_K"
        )
    return float(res["deff_H"])


def effective_sample_size(weights: np.ndarray) -> float:
    """Effective sample size = len(weights) / design_effect(weights).

    Kish (1965). Survey Sampling. Wiley.
    """
    return float(weights.shape[0] / design_effect(weights))
