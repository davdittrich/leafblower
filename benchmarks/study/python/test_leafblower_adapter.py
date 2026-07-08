"""benchmarks/study/python/test_leafblower_adapter.py -- WU-7 (leafblower-2ouc.8)
Contract-shape + golden pytest suite for
benchmarks/study/python/leafblower_adapter.py.

Usage (repo root, single-thread BLAS per CLAUDE.md determinism rule):
    cd python && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \\
        .venv/bin/python -m pytest ../benchmarks/study/python/test_leafblower_adapter.py -q
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import math  # noqa: E402
import sys  # noqa: E402
from pathlib import Path  # noqa: E402

import pytest  # noqa: E402

_THIS_DIR = Path(__file__).resolve().parent          # benchmarks/study/python
_STUDY_DIR = _THIS_DIR.parent                          # benchmarks/study
sys.path.insert(0, str(_THIS_DIR))
sys.path.insert(0, str(_STUDY_DIR / "common"))

from leafblower_adapter import (  # noqa: E402
    LBW_METHODS,
    LEAFBLOWER_PY_ADAPTERS,
    _REPO_ROOT,
    run_leafblower,
)
from metrics import KL_NATIVE_FAMILIES  # noqa: E402
from problem_io import load_problem_spec  # noqa: E402

BUILD = os.environ.get("LBW_BUILD_TAG", "native")

CONTRACT_KEYS = {
    "weights_ref", "iterations", "status", "converged", "error_message",
    "wall_time_s", "peak_rss_bytes",
}
STATUS_ENUM = {
    "converged", "no_conv", "infeasible", "bound_violation", "bad_arg",
    "budget", "stall", "error",
}


@pytest.fixture(scope="module")
def toy():
    return load_problem_spec(_STUDY_DIR / "spec" / "toy_inline.json")


def _read_weights(weights_ref: str):
    import pandas as pd
    return pd.read_parquet(_REPO_ROOT / weights_ref)["weight"].to_numpy()


# ---------------------------------------------------------------------------
# Contract-shape checks (all 9 methods)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("method", LBW_METHODS)
def test_contract_shape(toy, method):
    res = run_leafblower(toy, method, BUILD)
    assert set(res.keys()) == CONTRACT_KEYS, f"{method}: unexpected key set {set(res.keys())}"
    assert res["status"] in STATUS_ENUM, f"{method}: status {res['status']!r} not in harmonized enum"
    assert (_REPO_ROOT / res["weights_ref"]).exists(), f"{method}: weights_ref missing on disk"
    w = _read_weights(res["weights_ref"])
    assert len(w) == len(toy["data"]), f"{method}: weights_ref has wrong row count"
    assert math.isfinite(res["wall_time_s"]) and res["wall_time_s"] > 0
    assert res["peak_rss_bytes"] is None or res["peak_rss_bytes"] > 0
    assert isinstance(res["converged"], bool)


# ---------------------------------------------------------------------------
# LEAFBLOWER_PY_ADAPTERS registry: 9 entries, leafblower_<method>_py keys
# ---------------------------------------------------------------------------

def test_registry_has_nine_entries():
    assert len(LEAFBLOWER_PY_ADAPTERS) == 9
    assert set(LEAFBLOWER_PY_ADAPTERS.keys()) == {f"leafblower_{m}_py" for m in LBW_METHODS}


def test_registry_entry_is_callable_and_contract_shaped(toy):
    res = LEAFBLOWER_PY_ADAPTERS["leafblower_raking_py"](toy, BUILD)
    assert set(res.keys()) == CONTRACT_KEYS


# ---------------------------------------------------------------------------
# Golden: toy_inline hand-derived closed form.
#
# toy_inline (spec/toy_inline.json): single margin "grp", 2 disjoint groups
# A={row1,row2}, B={row3,row4}, design_weights=(1,1,2,2), targets A=B=0.5.
# leafblower enforces Sigma(w) = n = 4 (row count) at exit, regardless of
# sum(design_weights)=6.
#
# For a single margin split into disjoint groups, the constraint set is two
# independent linear equalities sum_A w_i = 0.5*4 = 2 and sum_B w_i =
# 0.5*4 = 2. Any KL/entropic-family projection from a starting point d onto
# a per-group-sum-only constraint is the UNIQUE multiplicative rescale
# within each group: w_i = d_i * (target_group_mass / d_group_sum). Here
# scale_A = 2 / (1+1) = 1, scale_B = 2 / (2+2) = 0.5, so
#   w = (1*1, 1*1, 2*0.5, 2*0.5) = (1, 1, 1, 1) EXACTLY.
# This closed form holds for every solver in KL_NATIVE_FAMILIES
# (metrics.py): oris, oris_soft, raking, sinkhorn, greenkhorn, newton_kl.
# logit/chebyshev/greg use a different divergence with no assumed common
# closed form here, so they get a looser (but still fully derived) check:
# convergence achieved, Sigma(w)=n, bounds respected.
# ---------------------------------------------------------------------------

KL_NATIVE_METHODS = [m for m in LBW_METHODS if m in KL_NATIVE_FAMILIES]
OTHER_METHODS = [m for m in LBW_METHODS if m not in KL_NATIVE_FAMILIES]


@pytest.mark.parametrize("method", KL_NATIVE_METHODS)
def test_golden_kl_native_exact_projection(toy, method):
    res = run_leafblower(toy, method, BUILD)
    assert res["converged"] is True, f"{method}: {res}"
    w = _read_weights(res["weights_ref"])
    assert w == pytest.approx([1.0, 1.0, 1.0, 1.0], abs=1e-4), f"{method}: weights {w}"


@pytest.mark.parametrize("method", OTHER_METHODS)
def test_golden_other_families_generic_sanity(toy, method):
    res = run_leafblower(toy, method, BUILD)
    assert res["converged"] is True, f"{method}: {res}"
    w = _read_weights(res["weights_ref"])
    assert w.sum() == pytest.approx(4.0, abs=1e-4)
    assert (w >= toy["bounds"]["min"] - 1e-9).all()
    assert (w <= toy["bounds"]["max"] + 1e-9).all()


# ---------------------------------------------------------------------------
# Error classification: bad_arg (ValueError path -- pre-validation error,
# max_weight < min_weight raised before the solver ever runs).
# ---------------------------------------------------------------------------

def test_error_classification_bad_arg(toy):
    bad_problem = dict(toy)
    bad_problem["bounds"] = {"min": 10, "max": 0}
    res = run_leafblower(bad_problem, "raking", BUILD)
    assert res["status"] == "bad_arg"
    assert res["converged"] is False
    assert (_REPO_ROOT / res["weights_ref"]).exists()
    w = _read_weights(res["weights_ref"])
    assert len(w) == len(toy["data"])
    import numpy as np
    assert np.all(np.isnan(w))
    assert res["error_message"]
