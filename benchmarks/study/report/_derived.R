# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# _derived.R — computes small per-figure aggregates from the raw parquets
# and writes them to report/_derived/*.parquet (gitignored, per invariant I3
# — derived artifacts are never committed, only regenerated).
#
# This is a THIN STUB: it wires the source() + output-dir plumbing and
# documents the pattern each downstream figure WU should follow. Add one
# aggregate per figure/table as needed — do NOT grow this file into a
# monolith; split per-figure aggregates into their own named objects below.
#
# Read-only on benchmarks/study/results/*.parquet — no metric recompute,
# only grouping/summarising of already-recomputed columns from _data.R's
# read_metrics()/read_agreement()/read_runs().

# Resolve this module's own directory regardless of caller's cwd (repo root
# or benchmarks/study/report/ — matches the resolution strategy in _data.R).
.candidate_report_dirs <- c(
  "benchmarks/study/report",
  "study/report",
  "."
)
.this_dir <- Find(function(d) file.exists(file.path(d, "_data.R")), .candidate_report_dirs)
if (is.null(.this_dir)) stop("_derived.R: could not locate report/_data.R from working directory: ", getwd())

source(file.path(.this_dir, "_data.R"))

DERIVED_DIR <- file.path(.this_dir, "_derived")
dir.create(DERIVED_DIR, showWarnings = FALSE, recursive = TRUE)

#' Write one derived aggregate to `_derived/<name>.parquet`.
#'
#' Pattern for downstream WUs adding a figure-specific aggregate:
#'   agg <- read_metrics() |> build_series() |> group_by(...) |> summarise(...)
#'   write_derived(agg, "fig_workprecision_kl")
write_derived <- function(df, name) {
  arrow::write_parquet(df, file.path(DERIVED_DIR, paste0(name, ".parquet")))
  invisible(df)
}

# --- example / smoke-test aggregate (kept minimal; downstream WUs add more) --
# One row per (solver, tier, build_role): median wall time on the headline
# build selection, demonstrating the write_derived() pattern end-to-end.
# Gated behind an env var so simply sourcing this file never has side effects
# (mirrors the "sources CLEAN, no error" DoD requirement for _data.R).
if (identical(Sys.getenv("LBW_REPORT_DERIVED_SMOKE"), "1")) {
  metrics <- read_metrics()
  agg <- metrics |>
    build_series(variant = "headline") |>
    dplyr::mutate(tier = tier_of(solver)) |>
    dplyr::group_by(solver, tier, build_role) |>
    dplyr::summarise(wall_time_s_median = median_na(wall_time_s_median), .groups = "drop")
  write_derived(agg, "_smoke_wall_time_by_solver")
}
