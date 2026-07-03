"""CR-E13 (leafblower-5ye4.13): a long-format target with a duplicated (variable,
level) row is ambiguous; _parse_target used dict(zip(...)) which silently kept only
the last. Both R and Python now raise on duplicate levels."""

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import warnings

import numpy as np
import pandas as pd
import pytest

from leafblower import harvest


def test_duplicate_level_target_rejected():
    n = 100
    df = pd.DataFrame({"a": np.random.RandomState(1).choice(["x", "y"], n)})
    # level 'x' appears twice for variable 'a'
    tdf = pd.DataFrame({"variable": ["a", "a", "a"],
                        "level": ["x", "y", "x"],
                        "proportion": [0.5, 0.5, 0.3]})
    with pytest.raises(ValueError, match="duplicate level"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            harvest(df, tdf, attach_weights=False)


def test_nan_level_duplicates_rejected():
    # NaN==NaN parity with R's anyDuplicated: two NaN-level rows are a duplicate.
    n = 100
    df = pd.DataFrame({"a": np.random.RandomState(1).choice(["x", "y"], n)})
    tdf = pd.DataFrame({"variable": ["a", "a"], "level": [np.nan, np.nan],
                        "proportion": [0.5, 0.5]})
    with pytest.raises(ValueError, match="duplicate level"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            harvest(df, tdf, attach_weights=False)


def test_unique_levels_still_accepted():
    n = 100
    df = pd.DataFrame({"a": np.random.RandomState(1).choice(["x", "y"], n)})
    tdf = pd.DataFrame({"variable": ["a", "a"], "level": ["x", "y"],
                        "proportion": [0.5, 0.5]})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        out = harvest(df, tdf, attach_weights=False)
    assert len(out["weights"]) == n
    assert out["result"]["status"] == 0
