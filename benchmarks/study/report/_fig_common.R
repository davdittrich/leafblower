# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# _fig_common.R — shared plumbing + constants for report/fig_work_precision.R
# and report/fig_walltime.R (WU 2ouc.13.2, RQ1 work-precision + wall-time
# figures). Sources _data.R (shared tier_of/stale_of/build_series/theme_lbw
# infra) and adds the small amount of encoding logic genuinely reused by
# BOTH figure scripts: tier color palette, build-variant shape palette, the
# log-axis display floor, the Pareto-envelope helper, and the stale-package
# label suffix.
#
# Read-only on results/*.parquet (via _data.R); no metric recompute here —
# these are display-layer helpers only (axis floors, orderings, palettes).

.candidate_report_dirs <- c(
  "benchmarks/study/report",
  "study/report",
  "."
)
.this_dir <- Find(function(d) file.exists(file.path(d, "_data.R")), .candidate_report_dirs)
if (is.null(.this_dir)) stop("_fig_common.R: could not locate report/_data.R from working directory: ", getwd())

source(file.path(.this_dir, "_data.R"))
suppressPackageStartupMessages(library(scales))

FIG_DIR <- file.path(.this_dir, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# --- palettes (Okabe-Ito colorblind-safe) --------------------------------

TIER_COLORS <- c(
  leafblower = "#D55E00",
  Tier1      = "#0072B2",
  Tier2      = "#009E73",
  Tier3      = "#CC79A7",
  other      = "grey50"
)

# build_role (from build_series(variant = "all")): "competitor" = single-
# build competitor run; "headline" = leafblower portable build (CRAN/PyPI
# default, the headline series); "delta" = leafblower native
# (-O3 -march=native) tuned-build overlay. Shapes pair headline/delta as
# filled-vs-hollow triangles so the same solver's two builds are visually
# linked; competitors get a filled circle.
BUILD_SHAPES <- c(competitor = 16, headline = 17, delta = 2)

# --- log-axis display floor ----------------------------------------------

#' Floor non-positive values for safe log10 display.
#'
#' Some native-divergence metrics contain exact-zero or floating-point-noise
#' negative values (a solver converging to the reference to within rounding
#' error) that `scale_*_log10()` cannot render. This clamps them to just
#' below the smallest genuine positive value observed in `x`, so degenerate
#' points are visually distinguishable (bunched at the axis floor) without
#' being confused with real small-but-nonzero errors. Pure display floor —
#' does NOT alter any stored metric column.
#'
#' @param x numeric vector (a single family's error column).
#' @return numeric vector, same length, with non-positive entries clamped.
floor_for_log <- function(x) {
  pos <- x[is.finite(x) & x > 0]
  if (length(pos) == 0) return(x)
  floor_val <- min(pos) / 10
  pmax(x, floor_val)
}

# --- Pareto envelope -------------------------------------------------------

#' Empirical Pareto envelope over (error, time) points: for each error
#' threshold, the minimum time achieved by any point with error at or below
#' that threshold. Computed by sorting on error then taking a running
#' minimum of time (never eyeballed) — a point is retained iff it sets a
#' new running-minimum time.
#'
#' @param df data.frame with numeric columns `err` and `wall_time_s_median`.
#' @return the Pareto-optimal subset of `df`, ordered by `err` ascending,
#'   with an added `envelope_time` column (== wall_time_s_median at those
#'   rows, kept for clarity in geom_step()).
pareto_front <- function(df) {
  df <- df[order(df$err), , drop = FALSE]
  running_min <- cummin(df$wall_time_s_median)
  keep <- c(TRUE, diff(running_min) < 0)
  out <- df[keep, , drop = FALSE]
  out$envelope_time <- running_min[keep]
  out
}

# --- stale-package display label -------------------------------------------

#' Append the "⚠ historical baseline" marker to solver labels for stale
#' packages (`stale_of()`, DESIGN.md §7).
stale_label <- function(solver) {
  ifelse(stale_of(solver), paste0(solver, " ⚠"), solver)
}

# --- solver factor ordering (tier, then median wall time) -----------------

#' Order a one-row-per-solver data.frame's `solver` column into a factor,
#' grouped by tier (leafblower, Tier1..3, other) then ascending median time.
#' Avoids a `forcats` dependency (keeps _data.R's dependency-light policy).
order_solver_factor <- function(df, tier_col = "tier", time_col = "wall_time_s_median") {
  ord <- order(factor(df[[tier_col]], levels = c("leafblower", "Tier1", "Tier2", "Tier3", "other")),
               df[[time_col]])
  df$solver <- factor(df$solver, levels = unique(rev(df$solver[ord])))
  df
}
