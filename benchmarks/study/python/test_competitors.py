"""benchmarks/study/python/test_competitors.py -- WU-6 (leafblower-2ouc.7).

Tests for the uniform-contract Python competitor adapters in
`competitors.py`. Strict separation: this file imports nothing from
leafblower's own src/, r_bridge.cpp, R/, python/leafblower/.

Run with (single-thread BLAS, venv python, per CLAUDE.md determinism rule):
  cd python && OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    .venv/bin/python -m pytest ../benchmarks/study/python/test_competitors.py -q
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_THIS_DIR = Path(__file__).resolve().parent
_STUDY_DIR = _THIS_DIR.parent
sys.path.insert(0, str(_THIS_DIR))
sys.path.insert(0, str(_STUDY_DIR / "common"))

import competitors  # noqa: E402

with open(_STUDY_DIR / "spec" / "status_enum.json") as _f:
    _STATUS_ENUM = set(json.load(_f)["$defs"]["StatusEnum"]["enum"])

_CONTRACT_KEYS = {
    "weights_ref", "iterations", "status", "converged",
    "error_message", "wall_time_s", "peak_rss_bytes",
}

_REPO_ROOT = _STUDY_DIR.parent.parent


# ---------------------------------------------------------------------------
# Shared problem builders
# ---------------------------------------------------------------------------

def _make_problem(problem_id: str, bounds_min: float, bounds_max: float) -> dict:
    """10-row, K=2 synthetic problem. Observed cell counts (design weight 1
    each): (A,X)=4, (A,Y)=1, (B,X)=1, (B,Y)=4 -- i.e. seed table
    T = [[4,1],[1,4]] (genuine interaction, odds ratio 16, NOT independent).
    Targets m1={A:0.7,B:0.3}, m2={X:0.6,Y:0.4} require real 2-D adjustment
    of BOTH margins (row sums 5/5 -> target 7/3; col sums 5/5 -> target
    6/4) -- this is the same non-degenerate configuration used to hand-
    verify the IPF <-> Sinkhorn(K=T, eps=1) equivalence (see
    test_ot_extraction_matches_ipf_golden below).
    """
    m1 = ["A"] * 4 + ["A"] + ["B"] + ["B"] * 4
    m2 = ["X"] * 4 + ["Y"] + ["X"] + ["Y"] * 4
    df = pd.DataFrame({"m1": m1, "m2": m2})
    return dict(
        id=problem_id,
        data=df,
        design_weights=np.ones(len(df), dtype=np.float64),
        margins=["m1", "m2"],
        targets={"m1": {"A": 0.7, "B": 0.3}, "m2": {"X": 0.6, "Y": 0.4}},
        bounds={"min": bounds_min, "max": bounds_max},
        tol=1e-3,
        objective_families=["kl"],
        K=2,
    )


def _unbounded_problem(problem_id: str) -> dict:
    return _make_problem(problem_id, 0.0, math.inf)


def _bounded_problem(problem_id: str) -> dict:
    # IPF-converged weights on this table lie in ~[0.54, 1.55] (see
    # test_ot_extraction_matches_ipf_golden) -- [0, 10] is finite but slack.
    return _make_problem(problem_id, 0.0, 10.0)


def _single_margin_problem(problem_id: str) -> dict:
    """K=1 problem -- used to exercise the OT adapters' K==2 bad_arg guard
    (registry.json K_max=2 for pot_sinkhorn/pot_greenkhorn/ott_jax_sinkhorn)."""
    df = pd.DataFrame({"m1": ["A"] * 5 + ["B"] * 5})
    return dict(
        id=problem_id,
        data=df,
        design_weights=np.ones(len(df), dtype=np.float64),
        margins=["m1"],
        targets={"m1": {"A": 0.5, "B": 0.5}},
        bounds={"min": 0.0, "max": math.inf},
        tol=1e-3,
        objective_families=["kl"],
        K=1,
    )


# ---------------------------------------------------------------------------
# Contract-shape assertion helper
# ---------------------------------------------------------------------------

def _assert_contract_shape(res: dict, n: int) -> None:
    assert set(res.keys()) == _CONTRACT_KEYS, res.keys()
    assert isinstance(res["weights_ref"], str) and res["weights_ref"]
    wpath = _REPO_ROOT / res["weights_ref"]
    assert wpath.exists(), f"weights_ref does not exist on disk: {wpath}"
    wdf = pd.read_parquet(wpath)
    assert len(wdf) == n
    assert res["iterations"] is None or isinstance(res["iterations"], int)
    assert res["status"] in _STATUS_ENUM, res["status"]
    assert isinstance(res["converged"], bool)
    assert res["error_message"] is None or isinstance(res["error_message"], str)
    assert isinstance(res["wall_time_s"], float) and res["wall_time_s"] > 0.0
    assert isinstance(res["peak_rss_bytes"], int) and res["peak_rss_bytes"] >= 0


def _read_weights(res: dict) -> np.ndarray:
    wpath = _REPO_ROOT / res["weights_ref"]
    return pd.read_parquet(wpath)["weight"].to_numpy(dtype=np.float64)


# ---------------------------------------------------------------------------
# Critical guard: OT 2-margin reduction must reproduce a hand-solved IPF
# golden (DESIGN.md Section 2 / ticket leafblower-2ouc.7 mandatory guard).
# ---------------------------------------------------------------------------

def _hand_ipf_multipliers() -> np.ndarray:
    """Reference IPF fixed point for T=[[4,1],[1,4]], target_row=[7,3],
    target_col=[6,4] (matches _make_problem's table). Returns the 2x2 cell
    multiplier matrix mult[i,j] = converged_ij / T_ij."""
    T = np.array([[4.0, 1.0], [1.0, 4.0]])
    target_row = np.array([7.0, 3.0])
    target_col = np.array([6.0, 4.0])
    M = T.copy()
    for _ in range(2000):
        M = M * (target_row / M.sum(axis=1))[:, None]
        M = M * (target_col / M.sum(axis=0))[None, :]
    return M / T


def _expected_ot_weights(problem: dict) -> np.ndarray:
    mult = _hand_ipf_multipliers()  # [row=m1 in {A,B}, col=m2 in {X,Y}]
    row_idx = {"A": 0, "B": 1}
    col_idx = {"X": 0, "Y": 1}
    d = problem["design_weights"]
    w = np.array([
        d[i] * mult[row_idx[r], col_idx[c]]
        for i, (r, c) in enumerate(zip(problem["data"]["m1"], problem["data"]["m2"]))
    ])
    return w * (d.sum() / w.sum())  # renormalize to sum(d), matches adapter convention


@pytest.mark.parametrize("solver_id", ["pot_sinkhorn", "pot_greenkhorn", "ott_jax_sinkhorn"])
def test_ot_extraction_matches_ipf_golden(solver_id):
    problem = _unbounded_problem(f"ot-ipf-golden-{solver_id}")
    res = competitors.ADAPTERS[solver_id](problem)
    _assert_contract_shape(res, len(problem["data"]))
    assert res["status"] == "converged", res["error_message"]
    assert res["converged"] is True

    expected = _expected_ot_weights(problem)
    got = _read_weights(res)
    np.testing.assert_allclose(got, expected, rtol=1e-3, atol=1e-6)


# ---------------------------------------------------------------------------
# Home-turf golden tests: one per adapter, each on its own canonical shape
# built from the same base problem (bounds adapted per adapter's family).
# ---------------------------------------------------------------------------

def _assert_sane_converged(res, problem):
    n = len(problem["data"])
    _assert_contract_shape(res, n)
    assert res["status"] == "converged", (res["status"], res["error_message"])
    assert res["converged"] is True, res["error_message"]
    w = _read_weights(res)
    assert np.all(np.isfinite(w))
    assert np.all(w > 0.0)
    np.testing.assert_allclose(w.sum(), problem["design_weights"].sum(), rtol=1e-6)


def test_ipfn_home_turf_golden():
    problem = _unbounded_problem("ipfn-golden")
    res = competitors.run_ipfn(problem)
    _assert_sane_converged(res, problem)


def test_weightipy_home_turf_golden():
    problem = _unbounded_problem("weightipy-golden")
    res = competitors.run_weightipy(problem)
    _assert_sane_converged(res, problem)


def test_svy_home_turf_golden():
    problem = _unbounded_problem("svy-golden")
    res = competitors.run_svy(problem)
    _assert_sane_converged(res, problem)


def test_balance_home_turf_golden():
    problem = _unbounded_problem("balance-golden")
    res = competitors.run_balance(problem)
    _assert_sane_converged(res, problem)


def test_scipy_trust_constr_home_turf_golden():
    problem = _bounded_problem("scipy-golden")
    res = competitors.run_scipy_trust_constr(problem)
    _assert_sane_converged(res, problem)


def test_cvxpy_linf_home_turf_golden():
    problem = _bounded_problem("cvxpy-golden")
    res = competitors.run_cvxpy_linf(problem)
    _assert_sane_converged(res, problem)


def test_pot_sinkhorn_home_turf_golden():
    problem = _unbounded_problem("pot-sinkhorn-golden")
    res = competitors.run_pot_sinkhorn(problem)
    _assert_sane_converged(res, problem)


def test_pot_greenkhorn_home_turf_golden():
    problem = _unbounded_problem("pot-greenkhorn-golden")
    res = competitors.run_pot_greenkhorn(problem)
    _assert_sane_converged(res, problem)


def test_ott_jax_sinkhorn_home_turf_golden():
    problem = _unbounded_problem("ott-jax-golden")
    res = competitors.run_ott_jax_sinkhorn(problem)
    _assert_sane_converged(res, problem)


# ---------------------------------------------------------------------------
# bad_arg guards
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("solver_id", ["ipfn", "weightipy", "svy", "balance"])
def test_bounded_problem_is_bad_arg_for_unbounded_only_adapters(solver_id):
    problem = _bounded_problem(f"bad-arg-bounds-{solver_id}")
    res = competitors.ADAPTERS[solver_id](problem)
    _assert_contract_shape(res, len(problem["data"]))
    assert res["status"] == "bad_arg"
    assert res["converged"] is False
    w = _read_weights(res)
    assert np.all(np.isnan(w))


@pytest.mark.parametrize("solver_id", ["pot_sinkhorn", "pot_greenkhorn", "ott_jax_sinkhorn"])
def test_k_not_two_is_bad_arg_for_ot_adapters(solver_id):
    problem = _single_margin_problem(f"bad-arg-k-{solver_id}")
    res = competitors.ADAPTERS[solver_id](problem)
    _assert_contract_shape(res, len(problem["data"]))
    assert res["status"] == "bad_arg"
    assert res["converged"] is False
