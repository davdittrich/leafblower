# STUDY-BRANCH-ONLY-DO-NOT-MERGE
"""benchmarks/study/analysis/aggregate_variance_cost.py -- WU (leafblower-2ouc.14.3)

Aggregates the per-replicate COLD re-solve shard parquets written by the
variance driver (leafblower-2ouc.14.2, R: variance_cost_run.R) into one row
per (solver, problem): B, replicate wall-time median/min/max, convergence
counts, total_variance_cost, and inflation_factor relative to the frozen
single-solve baseline in results/runs.parquet.

Reuses aggregate_metrics.py's registry helpers (load_registry,
family_primary_for) verbatim -- no new family-mapping logic here.

Analysis/read-only: does NOT edit any leafblower-package file, or any
benchmarks/study/common/ file, or aggregate_metrics.py itself.

Run (repo root, single-thread BLAS -- see /home/dd/Gemini/leafblower/CLAUDE.md):
    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        python/.venv/bin/python benchmarks/study/analysis/aggregate_variance_cost.py
"""

from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import sys  # noqa: E402
from pathlib import Path  # noqa: E402
from typing import Any  # noqa: E402

import polars as pl  # noqa: E402

_THIS_DIR = Path(__file__).resolve().parent          # benchmarks/study/analysis
STUDY_DIR = _THIS_DIR.parent                          # benchmarks/study
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

from aggregate_metrics import family_primary_for, load_registry  # noqa: E402

RESULTS_DIR = STUDY_DIR / "results"
RUNS_PATH = RESULTS_DIR / "runs.parquet"
SHARDS_DIR = RESULTS_DIR / "_variance_shards"
VARIANCE_COST_OUT = RESULTS_DIR / "variance_cost.parquet"


def aggregate_variance_cost(
    shards_df: pl.DataFrame, runs_df: pl.DataFrame, registry: dict[str, dict[str, Any]]
) -> pl.DataFrame:
    """One row per (solver, problem) group of shard replicate rows.

    inflation_factor = per_replicate_wall_median / base_solve_median, where
    base_solve_median is the median wall_time_s of runs_df's converged==True
    rows for that (solver, problem). NULL (never computed/divided) when the
    base cell is absent from runs_df OR present with zero converged rows --
    a base with no valid solve time has no meaningful ratio.
    """
    # base_solve_median per (solver, problem): median over CONVERGED base
    # reps only. Cells with zero converged base reps are simply absent from
    # this table (join below leaves inflation_factor null for them) -- this
    # is how we avoid ever dividing by a failed-base median.
    base_medians = (
        runs_df.filter(pl.col("converged"))
        .group_by(["solver", "problem"])
        .agg(pl.col("wall_time_s").median().alias("base_solve_median"))
    )

    cells = (
        shards_df.sort(["solver", "problem", "rep"])
        .group_by(["solver", "problem"], maintain_order=True)
        .agg(
            pl.len().alias("B"),
            pl.col("wall_time_s").median().alias("per_replicate_wall_median"),
            pl.col("wall_time_s").min().alias("per_replicate_wall_min"),
            pl.col("wall_time_s").max().alias("per_replicate_wall_max"),
            pl.col("converged").sum().alias("n_converged"),
        )
        .with_columns(
            (pl.col("B") - pl.col("n_converged")).alias("n_nonconv"),
            (pl.col("B") * pl.col("per_replicate_wall_median")).alias("total_variance_cost"),
        )
    )

    cells = cells.join(base_medians, on=["solver", "problem"], how="left").with_columns(
        (pl.col("per_replicate_wall_median") / pl.col("base_solve_median")).alias("inflation_factor")
    ).drop("base_solve_median")

    families = [
        family_primary_for(registry[solver]) if solver in registry else None
        for solver in cells["solver"].to_list()
    ]
    cells = cells.with_columns(pl.Series("family", families))

    return cells.select(
        "solver", "problem", "family", "B",
        "per_replicate_wall_median", "per_replicate_wall_min", "per_replicate_wall_max",
        "total_variance_cost", "n_converged", "n_nonconv", "inflation_factor",
    )


def main() -> None:
    registry = load_registry()
    shard_paths = sorted(SHARDS_DIR.glob("*.parquet"))
    if not shard_paths:
        raise FileNotFoundError(f"no shard parquet files found under {SHARDS_DIR}")
    shards_df = pl.concat(
        [pl.read_parquet(p) for p in shard_paths],
        how="diagonal_relaxed",
    )
    runs_df = pl.read_parquet(RUNS_PATH)

    out = aggregate_variance_cost(shards_df, runs_df, registry)
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out.write_parquet(VARIANCE_COST_OUT)

    print(f"n_variance_cost_rows={out.height}")
    n_null_inflation = out.filter(pl.col("inflation_factor").is_null()).height
    print(f"n_null_inflation_factor={n_null_inflation}")


if __name__ == "__main__":
    main()
