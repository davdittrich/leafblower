# STUDY-BRANCH-ONLY-DO-NOT-MERGE
"""benchmarks/study/analysis/test_aggregate_variance_cost.py -- WU (leafblower-2ouc.14.3)

pytest tests for aggregate_variance_cost.py's pure aggregation function.
Builds small in-memory polars fixtures (no real run needed) covering:
  - a (solver, problem) cell whose base is present + converged in
    runs.parquet -> finite inflation_factor.
  - a (solver, problem) cell whose base is ABSENT from runs.parquet ->
    inflation_factor NULL, no crash.
  - a (solver, problem) cell whose base rows are ALL non-converged ->
    inflation_factor NULL, no crash (never divide by a failed-base median).

Run (repo root, single-thread BLAS -- see CLAUDE.md):
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        python/.venv/bin/python -m pytest benchmarks/study/analysis/test_aggregate_variance_cost.py -q
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import math
import sys
from pathlib import Path

import polars as pl
import pytest

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

import aggregate_variance_cost as avc  # noqa: E402

# solver_a / problem_p: base present + converged in runs.parquet.
# solver_b / problem_p: base ABSENT from runs.parquet entirely.
# solver_c / problem_p: base present but ALL non-converged (all-failed).
SOLVER_A = "leafblower_raking_r"
SOLVER_B = "leafblower_greg_r"
SOLVER_C = "leafblower_logit_r"
PROBLEM = "toy_problem"


@pytest.fixture()
def registry() -> dict:
    entry = {"families": ["kl"], "method_class": "raking", "bounds": "both"}
    return {SOLVER_A: entry, SOLVER_B: entry, SOLVER_C: entry}


def _shard_rows(solver: str, problem: str, wall_times: list[float], converged: list[bool]) -> list[dict]:
    assert len(wall_times) == len(converged)
    return [
        dict(
            solver=solver, problem=problem, rep=i + 1, wall_time_s=wt,
            status=("converged" if conv else "error"), converged=conv,
            iterations=10, repwt_fingerprint="deadbeef",
        )
        for i, (wt, conv) in enumerate(zip(wall_times, converged))
    ]


@pytest.fixture()
def shards_df() -> pl.DataFrame:
    rows: list[dict] = []
    # solver_a: B=5, mixed converged (4 converged, 1 not), replicate wall times.
    rows += _shard_rows(
        SOLVER_A, PROBLEM,
        wall_times=[0.10, 0.11, 0.12, 0.13, 0.50],
        converged=[True, True, True, True, False],
    )
    # solver_b: B=5, all converged -- base is absent from runs.parquet.
    rows += _shard_rows(
        SOLVER_B, PROBLEM,
        wall_times=[0.20, 0.21, 0.22, 0.23, 0.24],
        converged=[True, True, True, True, True],
    )
    # solver_c: B=5, mixed converged -- base rows in runs.parquet are ALL non-converged.
    rows += _shard_rows(
        SOLVER_C, PROBLEM,
        wall_times=[0.30, 0.31, 0.32, 0.33, 0.34],
        converged=[True, True, False, True, True],
    )
    return pl.DataFrame(rows, infer_schema_length=None)


@pytest.fixture()
def runs_df() -> pl.DataFrame:
    rows: list[dict] = []
    # solver_a base: present, converged reps -- median wall_time_s = 1.0.
    for wt, conv in zip([0.9, 1.0, 1.1], [True, True, True]):
        rows.append(dict(
            solver=SOLVER_A, problem=PROBLEM, thread=1, build="portable", rep=1,
            wall_time_s=wt, status=("converged" if conv else "error"),
            converged=conv, weights_ref=None, peak_rss_bytes=1000, iterations=5,
        ))
    # solver_b base: intentionally ABSENT (no rows at all).
    # solver_c base: present but ALL non-converged.
    for wt in [2.0, 2.1, 2.2]:
        rows.append(dict(
            solver=SOLVER_C, problem=PROBLEM, thread=1, build="portable", rep=1,
            wall_time_s=wt, status="error", converged=False,
            weights_ref=None, peak_rss_bytes=1000, iterations=5,
        ))
    return pl.DataFrame(rows, infer_schema_length=None)


def test_one_row_per_solver_problem(shards_df, runs_df, registry):
    out = avc.aggregate_variance_cost(shards_df, runs_df, registry)
    pairs = set(zip(out["solver"].to_list(), out["problem"].to_list()))
    assert pairs == {(SOLVER_A, PROBLEM), (SOLVER_B, PROBLEM), (SOLVER_C, PROBLEM)}
    assert out.height == 3


def test_expected_columns_present(shards_df, runs_df, registry):
    out = avc.aggregate_variance_cost(shards_df, runs_df, registry)
    expected = {
        "solver", "problem", "family", "B",
        "per_replicate_wall_median", "per_replicate_wall_min", "per_replicate_wall_max",
        "total_variance_cost", "n_converged", "n_nonconv", "inflation_factor",
    }
    assert expected <= set(out.columns)


def test_family_sourced_from_registry(shards_df, runs_df, registry):
    out = avc.aggregate_variance_cost(shards_df, runs_df, registry)
    row = out.filter(pl.col("solver") == SOLVER_A).row(0, named=True)
    assert row["family"] == "kl"


def test_B_and_convergence_counts(shards_df, runs_df, registry):
    out = avc.aggregate_variance_cost(shards_df, runs_df, registry)
    row_a = out.filter(pl.col("solver") == SOLVER_A).row(0, named=True)
    assert row_a["B"] == 5
    assert row_a["n_converged"] == 4
    assert row_a["n_nonconv"] == 1
    assert row_a["B"] == row_a["n_converged"] + row_a["n_nonconv"]


def test_total_variance_cost_equals_B_times_median(shards_df, runs_df, registry):
    out = avc.aggregate_variance_cost(shards_df, runs_df, registry)
    row_a = out.filter(pl.col("solver") == SOLVER_A).row(0, named=True)
    expected_median = sorted([0.10, 0.11, 0.12, 0.13, 0.50])[2]
    assert row_a["per_replicate_wall_median"] == pytest.approx(expected_median)
    assert row_a["total_variance_cost"] == pytest.approx(row_a["B"] * expected_median)


def test_inflation_factor_finite_when_base_present_and_converged(shards_df, runs_df, registry):
    out = avc.aggregate_variance_cost(shards_df, runs_df, registry)
    row_a = out.filter(pl.col("solver") == SOLVER_A).row(0, named=True)
    base_median = 1.0  # median([0.9, 1.0, 1.1])
    expected = row_a["per_replicate_wall_median"] / base_median
    assert row_a["inflation_factor"] is not None
    assert math.isfinite(row_a["inflation_factor"])
    assert row_a["inflation_factor"] == pytest.approx(expected)


def test_inflation_factor_null_when_base_absent(shards_df, runs_df, registry):
    out = avc.aggregate_variance_cost(shards_df, runs_df, registry)
    row_b = out.filter(pl.col("solver") == SOLVER_B).row(0, named=True)
    assert row_b["inflation_factor"] is None


def test_inflation_factor_null_when_base_all_failed(shards_df, runs_df, registry):
    out = avc.aggregate_variance_cost(shards_df, runs_df, registry)
    row_c = out.filter(pl.col("solver") == SOLVER_C).row(0, named=True)
    assert row_c["inflation_factor"] is None
    # sanity: base cell for solver_c IS present in runs_df, just all-failed.
    base_c = runs_df.filter(
        (pl.col("solver") == SOLVER_C) & (pl.col("problem") == PROBLEM)
    )
    assert base_c.height > 0
    assert base_c.filter(pl.col("converged")).height == 0


def test_no_exception_raised_end_to_end(shards_df, runs_df, registry):
    # Full pipeline (all three cases mixed) must not raise.
    out = avc.aggregate_variance_cost(shards_df, runs_df, registry)
    assert out.height == 3
