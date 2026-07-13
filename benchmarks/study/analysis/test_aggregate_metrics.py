# STUDY-BRANCH-ONLY-DO-NOT-MERGE
"""benchmarks/study/analysis/test_aggregate_metrics.py -- WU-12a smoke test.

Script-style (run directly, not via pytest collection) -- asserts that
aggregate_metrics.py's output reuses common/metrics.py's primitives
correctly. Run AFTER aggregate_metrics.py has produced
results/metrics.parquet + results/agreement.parquet:

    OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
        python/.venv/bin/python benchmarks/study/analysis/test_aggregate_metrics.py

Checks:
  (a) one competitor cell's DEFF_kish/marg_kl_mean/marg_kl_max == a FRESH
      direct metrics.compute_metrics() call on the same weight vector
      (spies the reuse -- proves aggregate_metrics.py doesn't reimplement
      the formula, just calls compute_metrics with the same inputs).
  (b) >=1 strictly-convex-family R<->Py agreement pair exists with
      pearson in [-1, 1].
  (c) >=1 minimax agreement pair in mode=='objective_value' with
      obj_val_a == that cell's margin_linf (from metrics.parquet).
  (d) native_div is finite for a chi2 cell AND a logit cell with finite
      bounds.
  (e) metrics.parquet is non-empty and has the expected columns.
"""

from __future__ import annotations

import math
import os
import sys
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import polars as pl  # noqa: E402

_THIS_DIR = Path(__file__).resolve().parent
if str(_THIS_DIR) not in sys.path:
    sys.path.insert(0, str(_THIS_DIR))

import aggregate_metrics as am  # noqa: E402
import metrics  # noqa: E402  (benchmarks/study/common, on sys.path via aggregate_metrics)

EXPECTED_METRICS_COLUMNS = {
    "solver", "problem", "thread", "build", "status", "rep_count",
    "wall_time_s_median", "wall_time_s_min", "wall_time_s_max",
    "peak_rss_bytes", "iterations", "method_class", "family",
    "native_div", "native_div_kind",
    "n", "W", "marg_kl_mean", "marg_kl_max", "marg_kl_divergent",
    "margin_linf", "margin_l1", "ESS", "DEFF_kish", "DEFF_g", "g_weighted",
    "weight_kl", "wmin", "wmed", "wmax",
    "bound_viol_count", "bound_viol_max", "bound_viol_mean",
}


def check_a(metrics_df: pl.DataFrame, runs_df: pl.DataFrame) -> None:
    """(a) competitor DEFF_kish/marg_kl == a fresh direct compute_metrics call."""
    competitor_rows = metrics_df.filter(~pl.col("solver").str.starts_with("leafblower"))
    assert competitor_rows.height > 0, "no competitor rows in metrics.parquet"
    cell = competitor_rows.row(0, named=True)

    run_row = (
        runs_df.filter(
            (pl.col("solver") == cell["solver"]) & (pl.col("problem") == cell["problem"])
            & (pl.col("thread") == cell["thread"]) & (pl.col("build") == cell["build"])
        )
        .row(0, named=True)
    )
    w = am.load_real_weight(run_row["weights_ref"])
    assert w is not None, f"could not reload real weights for {cell['solver']}/{cell['problem']}"

    problem = am.load_problem(cell["problem"])
    bnds = am.bounds_for_metrics(problem["bounds"])
    fresh = metrics.compute_metrics(
        w, problem["groups"], problem["targets"], d=problem["design_weights"],
        bounds=bnds, family=cell["family"],
    )
    assert fresh["DEFF_kish"] == cell["DEFF_kish"], (
        f"DEFF_kish mismatch: fresh={fresh['DEFF_kish']} stored={cell['DEFF_kish']}"
    )
    assert fresh["marg_kl_mean"] == cell["marg_kl_mean"], (
        f"marg_kl_mean mismatch: fresh={fresh['marg_kl_mean']} stored={cell['marg_kl_mean']}"
    )
    assert fresh["marg_kl_max"] == cell["marg_kl_max"], (
        f"marg_kl_max mismatch: fresh={fresh['marg_kl_max']} stored={cell['marg_kl_max']}"
    )
    print(f"(a) PASS -- {cell['solver']}/{cell['problem']} DEFF_kish={cell['DEFF_kish']!r} "
          f"marg_kl_mean={cell['marg_kl_mean']!r} bit-identical to a fresh compute_metrics() call")


def check_b(agreement_df: pl.DataFrame) -> None:
    """(b) >=1 strictly-convex R<->Py pair with pearson in [-1, 1]."""
    convex = agreement_df.filter(
        (pl.col("mode") == "weight_vector") & (pl.col("family").is_in(list(am.STRICTLY_CONVEX)))
    )
    rp_pairs = convex.filter(
        pl.struct(["solver_a", "solver_b"]).map_elements(
            lambda s: am.base_solver_name(s["solver_a"]) == am.base_solver_name(s["solver_b"])
            and s["solver_a"] != s["solver_b"],
            return_dtype=pl.Boolean,
        )
    )
    assert rp_pairs.height > 0, "no strictly-convex R<->Py agreement pair found"
    # A constant weight vector (zero variance) yields NaN pearson (corrcoef /0);
    # require >=1 pair with a genuinely finite pearson in [-1, 1].
    finite = rp_pairs.filter(
        pl.col("pearson").is_finite() & (pl.col("pearson") >= -1.0) & (pl.col("pearson") <= 1.0)
    )
    assert finite.height > 0, "no strictly-convex R<->Py pair with finite pearson in [-1, 1]"
    row = finite.row(0, named=True)
    print(f"(b) PASS -- {row['solver_a']}({row['build_a']}) <-> {row['solver_b']}({row['build_b']}) "
          f"on {row['problem']} ({row['family']}): pearson={row['pearson']!r} "
          f"({rp_pairs.height} R<->Py convex pairs, {finite.height} with finite pearson)")


def check_c(agreement_df: pl.DataFrame, metrics_df: pl.DataFrame) -> None:
    """(c) >=1 minimax pair in objective_value mode with obj_val==cell's margin_linf."""
    minimax_pairs = agreement_df.filter(
        (pl.col("family") == "minimax") & (pl.col("mode") == "objective_value")
    )
    assert minimax_pairs.height > 0, "no minimax objective_value agreement pair found"
    row = minimax_pairs.row(0, named=True)
    cell = metrics_df.filter(
        (pl.col("solver") == row["solver_a"]) & (pl.col("build") == row["build_a"])
        & (pl.col("problem") == row["problem"])
    ).row(0, named=True)
    assert row["obj_val_a"] == cell["margin_linf"], (
        f"obj_val_a ({row['obj_val_a']}) != cell margin_linf ({cell['margin_linf']})"
    )
    print(f"(c) PASS -- {row['solver_a']}({row['build_a']}) vs {row['solver_b']}({row['build_b']}) "
          f"on {row['problem']}: obj_val_a={row['obj_val_a']!r} == margin_linf")


def check_d(metrics_df: pl.DataFrame) -> None:
    """(d) native_div finite for a chi2 cell AND a logit cell with finite bounds."""
    chi2_rows = metrics_df.filter(
        (pl.col("family") == "chi2") & pl.col("native_div").is_not_null()
    )
    assert chi2_rows.height > 0, "no chi2 cell with non-null native_div"
    chi2_cell = chi2_rows.row(0, named=True)
    assert math.isfinite(chi2_cell["native_div"]), f"chi2 native_div not finite: {chi2_cell['native_div']}"

    logit_rows = metrics_df.filter(
        (pl.col("family") == "logit") & pl.col("native_div").is_not_null()
    )
    assert logit_rows.height > 0, "no logit cell with finite bounds (non-null native_div)"
    logit_cell = logit_rows.row(0, named=True)
    assert math.isfinite(logit_cell["native_div"]), f"logit native_div not finite: {logit_cell['native_div']}"
    print(f"(d) PASS -- chi2 native_div={chi2_cell['native_div']!r} ({chi2_cell['solver']}/{chi2_cell['problem']}), "
          f"logit native_div={logit_cell['native_div']!r} ({logit_cell['solver']}/{logit_cell['problem']})")


def check_e(metrics_df: pl.DataFrame) -> None:
    """(e) metrics.parquet non-empty + has expected columns."""
    assert metrics_df.height > 0, "metrics.parquet is empty"
    missing = EXPECTED_METRICS_COLUMNS - set(metrics_df.columns)
    assert not missing, f"metrics.parquet missing expected columns: {missing}"
    print(f"(e) PASS -- metrics.parquet has {metrics_df.height} rows, all expected columns present")


def main() -> None:
    assert am.METRICS_OUT.exists(), f"{am.METRICS_OUT} not found -- run aggregate_metrics.py first"
    assert am.AGREEMENT_OUT.exists(), f"{am.AGREEMENT_OUT} not found -- run aggregate_metrics.py first"
    metrics_df = pl.read_parquet(am.METRICS_OUT)
    agreement_df = pl.read_parquet(am.AGREEMENT_OUT)
    runs_df = pl.read_parquet(am.RUNS_PATH)

    check_e(metrics_df)
    check_a(metrics_df, runs_df)
    check_b(agreement_df)
    check_c(agreement_df, metrics_df)
    check_d(metrics_df)

    print("ALL CHECKS PASS")


if __name__ == "__main__":
    try:
        main()
    except AssertionError as e:
        print(f"FAIL: {e}")
        sys.exit(1)
