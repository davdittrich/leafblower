# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# tbl_quality.R — RQ4 quality tables (median+quantile ESS + DEFF, per
# (family, solver)) + RQ5 minimax objective-value agreement TABLE (rel_gap,
# NEVER weight-vector correlation — DESIGN.md §6 "Blocker G": the L-inf LP
# optimum lies on a face, not a unique vertex, so two correct solvers return
# different weight vectors at identical achieved L-inf; minimax is judged on
# achieved-objective agreement only).
#
# Writes report/tables/quality.csv + report/tables/minimax_objval.csv.
# Read-only on results/metrics.parquet + results/agreement.parquet via
# _data.R — no metric recompute, only grouping/summarising of already-
# recomputed columns (ESS, DEFF_kish, DEFF_g, weight_kl, rel_gap all exist
# upstream; median/quantile here are display-layer aggregates across problem
# instances, not a metric recompute).

.candidate_report_dirs <- c("benchmarks/study/report", "study/report", ".")
.this_dir <- Find(function(d) file.exists(file.path(d, "_data.R")), .candidate_report_dirs)
if (is.null(.this_dir)) stop("tbl_quality.R: could not locate report/_data.R from working directory: ", getwd())
source(file.path(.this_dir, "_data.R"))
suppressPackageStartupMessages(library(dplyr))

TBL_DIR <- file.path(.this_dir, "tables")
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

q25 <- function(x) stats::quantile(x, probs = 0.25, na.rm = TRUE, names = FALSE)
q75 <- function(x) stats::quantile(x, probs = 0.75, na.rm = TRUE, names = FALSE)

#' Write `df` to TBL_DIR/<name>.csv and assert it landed non-empty.
write_table <- function(df, name) {
  path <- file.path(TBL_DIR, paste0(name, ".csv"))
  utils::write.csv(df, path, row.names = FALSE)
  stopifnot(
    "table was not written" = file.exists(path),
    "table is empty" = file.info(path)$size > 0
  )
  message("wrote ", path)
  invisible(df)
}

# --- RQ4 quality: median + quantile ESS / DEFF, per (family, solver) ------
#
# build_series(variant = "headline") is the DESIGN.md §5 chokepoint (one row
# per (solver, problem): leafblower restricted to `portable`, competitors'
# single build pass through) — avoids double-counting a leafblower solver's
# portable+native rows against single-build competitors when computing the
# per-solver median/quantile below.
metrics_hl <- read_metrics() |> build_series(variant = "headline")

stopifnot(
  "tbl_quality: metrics_hl not one row per (solver,problem) after headline build_series (thread>1?)" =
    !any(duplicated(metrics_hl[c("solver", "problem")]))
)

# DEFF_g (`1+CV^2(g)`, g_i = w_i/d_i) is reported ALONGSIDE DEFF_kish
# (`1+CV^2(w)`, the Kish weighting-loss = n/ESS — explicitly NOT the true
# design effect), not as a replacement: on d_i==1 problems the two coincide,
# but on d_i!=1 problems (g_weighted flags these; verified against
# metrics.parquet 2026-07-13: canonical_survey_apistrat + toy_inline in the
# current results) raw-w Kish DEFF conflates intentional base-design
# variance with calibration-injected variance (DESIGN.md Gap D). DEFF_g_med
# is computed over ONLY the g_weighted==TRUE rows of each (family, solver)
# group (median_na(numeric(0)) == NA when a solver has none — i.e. every
# problem it ran had d_i==1 — never a spurious 0/1).
quality <- metrics_hl |>
  mutate(tier = tier_of(solver), stale = stale_of(solver)) |>
  group_by(family, solver, tier, stale) |>
  summarise(
    n_problems       = dplyr::n(),
    ESS_median       = median_na(ESS),
    ESS_q25          = q25(ESS),
    ESS_q75          = q75(ESS),
    DEFF_kish_median = median_na(DEFF_kish),
    DEFF_kish_q25    = q25(DEFF_kish),
    DEFF_kish_q75    = q75(DEFF_kish),
    n_g_weighted     = sum(g_weighted),
    DEFF_g_median    = median_na(DEFF_g[g_weighted]),
    # weight_kl is this family's OWN objective for kl-family solvers (native)
    # and a non-native cross-family closeness-to-design diagnostic for every
    # other family (DESIGN.md §6 Gap B) — weight_kl_note carries the exact
    # caveat text verbatim from metrics.parquet, not re-derived here.
    weight_kl_median = median_na(weight_kl),
    weight_kl_note   = dplyr::first(weight_kl_note),
    .groups = "drop"
  ) |>
  arrange(family, solver)

stopifnot("quality table produced zero rows" = nrow(quality) > 0)
write_table(quality, "quality")

# --- RQ5 minimax: objective-value agreement TABLE (rel_gap) --------------
#
# mode == "objective_value" rows only (family == "minimax"); pearson/
# spearman/max_abs_diff/cosine are NA for these by construction
# (analysis/aggregate_metrics.py:276) and are never read here — this table
# is rel_gap only, never a vector-correlation on minimax weights.
agreement <- read_agreement()
.is_headline_build <- function(build_col) build_col %in% c("na", "portable")

minimax_filtered <- agreement |>
  filter(mode == "objective_value", family == "minimax") |>
  filter(.is_headline_build(build_a), .is_headline_build(build_b))

stopifnot(
  "tbl_quality minimax: >1 row per (solver_a,solver_b,problem) (thread>1?)" =
    !any(duplicated(minimax_filtered[c("solver_a", "solver_b", "problem")]))
)

minimax_objval <- minimax_filtered |>
  group_by(solver_a, solver_b) |>
  summarise(
    n_problems     = dplyr::n(),
    rel_gap_median = median_na(rel_gap),
    rel_gap_q25    = q25(rel_gap),
    rel_gap_q75    = q75(rel_gap),
    .groups = "drop"
  ) |>
  arrange(solver_a, solver_b)

stopifnot("minimax_objval table produced zero rows" = nrow(minimax_objval) > 0)
write_table(minimax_objval, "minimax_objval")
