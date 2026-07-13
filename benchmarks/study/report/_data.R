# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# _data.R — shared data-load + tier-label + build-selection + maintenance-
# status module for benchmarks/study/report/*.qmd.
#
# STRICT SEPARATION: read-only on benchmarks/study/results/*.parquet — no
# metric recompute. Sourced by report.qmd and _derived.R; keep it dependency-
# light (arrow, dplyr, ggplot2, scales only).

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(ggplot2)
})

# --- paths ------------------------------------------------------------
# Fallback: allow sourcing from repo root or from benchmarks/study/report/.
.candidate_results_dirs <- c(
  "benchmarks/study/results",
  "study/results",
  "../results",
  file.path(dirname(normalizePath(".")), "results")
)
.resolve_results_dir <- function() {
  for (d in .candidate_results_dirs) {
    if (dir.exists(d) && file.exists(file.path(d, "metrics.parquet"))) return(normalizePath(d))
  }
  stop("_data.R: could not locate benchmarks/study/results/ (metrics.parquet not found) ",
       "from working directory: ", getwd())
}
RESULTS_DIR <- .resolve_results_dir()

# --- readers ------------------------------------------------------------

#' Read the recomputed per-(solver,problem,thread,build) quality/timing metrics.
read_metrics <- function(path = file.path(RESULTS_DIR, "metrics.parquet")) {
  arrow::read_parquet(path)
}

#' Read the pairwise weight-agreement table (RQ5).
read_agreement <- function(path = file.path(RESULTS_DIR, "agreement.parquet")) {
  arrow::read_parquet(path)
}

#' Read the raw per-run harness log (timing reps, status, etc).
read_runs <- function(path = file.path(RESULTS_DIR, "runs.parquet")) {
  arrow::read_parquet(path)
}

# --- objective-family strata (DESIGN.md §6) -------------------------
strata <- c("kl", "chi2", "logit", "minimax")

# --- tier taxonomy (DESIGN.md §2 "Competitor discipline tiers", WU-12) ----
# Exact solver ids verified against benchmarks/study/registry.json
# (schema_version 2.0.0, 47 solvers total: 19 Tier1 + 7 Tier2 + 3 Tier3 + 18
# leafblower_*_{py,r}) on 2026-07-13.

.tier1_solvers <- c(
  "survey_calibrate_raking", "survey_calibrate_linear", "survey_calibrate_logit",
  "sampling_calib_linear", "sampling_calib_logit",
  "anesrake", "ipfr", "ipfn", "svy", "weightipy",
  "balance", "samplics", "regenesees", "icarus", "laeken",
  "gecal", "jointcalib", "nonprobsvy", "ebal"
)

.tier2_solvers <- c(
  "sbw", "weightit", "optweight_linf", "optweight_entropy",
  "pot_sinkhorn", "pot_greenkhorn", "ott_jax_sinkhorn"
)

.tier3_solvers <- c(
  "cvxpy_linf", "cvxr_reference", "scipy_trust_constr"
)

#' Map a solver id (character vector) to its competitor-discipline tier.
#'
#' @param solver character vector of registry solver ids (e.g. `metrics$solver`).
#' @return character vector, one of "leafblower", "Tier1", "Tier2", "Tier3",
#'   "other" (unrecognised id — never errors).
tier_of <- function(solver) {
  out <- rep("other", length(solver))
  out[grepl("^leafblower_", solver)] <- "leafblower"
  out[solver %in% .tier1_solvers] <- "Tier1"
  out[solver %in% .tier2_solvers] <- "Tier2"
  out[solver %in% .tier3_solvers] <- "Tier3"
  out
}

# --- maintenance status (DESIGN.md §7 "⚠ packages labelled 'historical baseline'") --
.stale_solvers <- c("anesrake", "ipfn", "regenesees")

#' Is `solver` in the ⚠ "historical baseline" / stale-maintenance set?
#' @return logical vector, same length as `solver`.
stale_of <- function(solver) {
  solver %in% .stale_solvers
}

# --- build-variant selection (DESIGN.md §5 "Build-flag parity") -----------
#
# Contract: leafblower ships two build variants (`build` ∈ {"portable",
# "native"}); competitors run once (`build == "na"`). A figure that pools
# ALL rows of `metrics`/`agreement` for a leafblower solver would silently
# double-count it against every single-build competitor. `build_series()`
# is the ONE chokepoint every WU MUST call before building a headline or
# native-delta figure/table — never filter on `build` ad hoc downstream.
#
# variant = "headline" (default): one row per solver — leafblower rows
#   restricted to build == "portable" (CRAN/PyPI-default optimisation,
#   matches competitors' shipped builds per DESIGN.md §5), competitor rows
#   (build == "na") pass through unchanged. This is what every RQ1-style
#   headline figure must use.
# variant = "native_delta": leafblower rows restricted to build == "native"
#   only (the "-O3 -march=native" tuned-build delta series); competitor
#   rows are dropped (there is no native competitor variant to compare).
# variant = "all": no filtering, but adds `build_role` so a caller can see
#   the composition (`"headline"` / `"delta"` / `"competitor"`) before
#   deciding how to facet/subset it further.
build_series <- function(df, variant = c("headline", "native_delta", "all")) {
  variant <- match.arg(variant)
  stopifnot("build" %in% names(df))

  df <- df |>
    mutate(build_role = case_when(
      build == "na"       ~ "competitor",
      build == "portable" ~ "headline",
      build == "native"   ~ "delta",
      TRUE ~ "other"
    ))

  switch(variant,
    headline     = df |> filter(build_role %in% c("competitor", "headline")),
    native_delta = df |> filter(build_role == "delta"),
    all          = df
  )
}

# --- small numeric helpers ------------------------------------------------

#' NA-safe median.
median_na <- function(x) stats::median(x, na.rm = TRUE)

#' Group-wise median of `value_col`, grouped by `...` (tidy-eval group cols).
grouped_median <- function(df, value_col, ...) {
  df |>
    group_by(...) |>
    summarise(median = median_na(.data[[value_col]]), .groups = "drop")
}

# --- shared ggplot theme (facet-based composition only; no patchwork) ----
theme_lbw <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.background = ggplot2::element_rect(fill = "grey92", colour = NA),
      strip.text = ggplot2::element_text(face = "bold"),
      legend.position = "bottom",
      plot.title.position = "plot"
    )
}
