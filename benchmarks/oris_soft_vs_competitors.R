# Honest Performance Gate — oris_soft vs. doc-named competitors
# Phase: 03-honest-performance-gate, Plan 01 (tracer slice)
#
# D-04: oris_soft is the single authoritative solver for the headline claim.
# D-07: the competitor set is drawn ONLY from docs/methods/oris.md's
#   "Practitioner implementations & use cases" table — this plan measures
#   `survey::calibrate` (raking), the remaining competitors arrive in plan 02.
# D-08: every number below is measured FRESH in this run. This script does
#   NOT read benchmarks/study/report/tables/*.csv or any other pre-aggregated
#   result file.
# D-05: no comparison against, mention of, or variable named after the
#   unreleased package the old "Nx faster than" framing used as its baseline.
#   That framing is exactly what this phase exists to retire — do not
#   reintroduce it, not even as a commented-out arm.
#
# CSV schema (frozen; plans 02-04 extend ROWS into it, not columns):
#   input_class, n, n_margins, n_categories, m_cell, m_cell_over_n,
#   max_weight, arm, wall_s, max_error, max_w, min_w, deff, n_eff,
#   iterations, ok, note

suppressPackageStartupMessages({
  library(leafblower)
  library(survey)
  library(bench)
})

# --- Determinism guard (CLAUDE.md single-thread BLAS protocol, ENFORCED) ---
# This is the first check the script performs, before any fixture or solve.
# Deliberately does NOT call Sys.setenv() to satisfy itself — the guard only
# has teeth if the caller's environment is what is checked.
require_single_thread_blas <- function() {
  vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
  vals <- Sys.getenv(vars, unset = NA_character_)
  bad  <- vars[is.na(vals) | vals != "1"]
  if (length(bad) > 0L) {
    stop(sprintf(
      "refusing to measure: %s must be set to \"1\" for a reproducible single-thread BLAS run (CLAUDE.md protocol); got {%s}",
      paste(bad, collapse = ", "),
      paste(sprintf("%s=%s", bad, ifelse(is.na(vals[bad]), "<unset>", vals[bad])), collapse = ", ")
    ), call. = FALSE)
  }
}

# Shared, arm-independent accuracy metric — applied identically to every
# arm's returned weight vector. Never reads a package's self-reported
# convergence number. Accumulates the achieved proportion directly (no
# difference-of-large-sums cancellation).
margin_max_error <- function(w, df, tgt) {
  Z <- sum(w)
  max(vapply(names(tgt), function(k) {
    T_k <- tgt[[k]]
    max(vapply(names(T_k), function(j) {
      abs(sum(w[df[[k]] == j]) / Z - T_k[[j]])
    }, numeric(1)))
  }, numeric(1)))
}

# Bound compliance is compared on one scale: sum(w) == n (leafblower's own
# exit convention, lbw::finalize_weights: a single pre-bounds scale to n,
# THEN bounds_mode dispatch). leafblower's returned weights already sit at
# this scale to double-precision (finalize_weights enforces it internally)
# and CLAUDE.md explicitly FORBIDS renormalizing AFTER unit-mode water-fill
# ("silently breaks the bounds_mode='unit' clamps") — confirmed empirically:
# re-scaling leafblower's own output pushed a weight sitting exactly at the
# max_weight bound to 3.0000000036, failing the max_w <= max_weight + 1e-10
# check this script itself asserts. So this helper is applied ONLY to the
# competitor arm below, whose weights (survey's epsilon-tolerance raking)
# do NOT already sit at sum(w) == n, to bring it onto the same scale.
normalize_to_n <- function(w, n) w * n / sum(w)

# One CSV row for one arm. Metrics are passed in pre-computed so each
# call site to margin_max_error() stays visible at the call site, not
# hidden behind this row-shaping helper.
arm_row <- function(input_class, n, n_margins, n_categories, m_cell, max_weight,
                     arm, wall_s, max_error, max_w, min_w, deff, n_eff,
                     iterations, ok, note) {
  data.frame(
    input_class = input_class, n = n, n_margins = n_margins,
    n_categories = n_categories, m_cell = m_cell, m_cell_over_n = m_cell / n,
    max_weight = max_weight, arm = arm, wall_s = wall_s, max_error = max_error,
    max_w = max_w, min_w = min_w, deff = deff, n_eff = n_eff,
    iterations = iterations, ok = ok, note = note,
    stringsAsFactors = FALSE
  )
}

# Runs both arms on one already-built input class (df + tgt) and returns
# the two CSV rows. Plan 02 adds further input classes and competitors by
# calling this same function (or its per-arm blocks) on new fixtures —
# the fixture/measurement split is what makes that extension additive.
run_input_class <- function(input_class, df, tgt, max_weight, n_categories) {
  n <- nrow(df)
  margin_cols <- names(tgt)
  n_margins   <- length(margin_cols)
  # A property of the fixture, identical for every arm — computed once.
  m_cell <- nrow(unique(df[margin_cols]))

  cat(sprintf("\n=== %s n=%d K=%d nj=%d m_cell=%d (m_cell/n=%.4f) ===\n",
              input_class, n, n_margins, n_categories, m_cell, m_cell / n))

  # --- Arm 1: leafblower_oris_soft ---
  # bounds_mode="unit" is load-bearing, not incidental: the default "cell"
  # enforces bounds only in cell aggregate, whereas survey/icarus/ReGenesees
  # all enforce a per-observation bound on the ratio g = w/d. With unit
  # design weights, bounds_mode="unit" + max_weight=3 corresponds exactly to
  # the competitors' bounds=c(0,3) on g — comparing a cell-aggregate
  # leafblower arm against per-observation competitors would be a different
  # optimisation problem and the timing comparison would be meaningless.
  # convergence is oris_soft's own CANONICAL default (R/harvest.R roxygen,
  # `convergence` param): metric="marginal_kl", rule="improvement", tol=0.001.
  # Passed explicitly (not convergence=list()) so the intent is visible at
  # the call site. Originally this arm passed the `absolute` shorthand,
  # which forces metric="max_err"/rule="threshold" — survey::calibrate's
  # own native stopping rule, not oris_soft's. That silently made oris_soft
  # stop the instant it first crossed the tolerance instead of continuing
  # to refine on its own plateau-detection rule, exactly the "stopped
  # early" confound this phase's determinism protocol exists to eliminate.
  # Fixed per tracer-checkpoint review; see 03-01-SUMMARY.md.
  lb_call <- function() {
    harvest(df, tgt, method = "oris_soft", max_weight = max_weight,
            bounds_mode = "unit", attach_weights = FALSE,
            convergence = list(metric = "marginal_kl", rule = "improvement", tol = 0.001))
  }
  bm_lb  <- bench::mark(run = lb_call(), iterations = 2, check = FALSE,
                         memory = FALSE, filter_gc = FALSE)
  w_lb   <- lb_call()
  res_lb <- attr(w_lb, "result")
  # NOT normalized — see the comment on normalize_to_n() above.
  w_lb_n <- as.numeric(w_lb)
  max_error_lb <- margin_max_error(w_lb_n, df, tgt)
  # status 0 = converged; 5 = plateau at constrained optimum (weights valid,
  # the bound is legitimately active) — both are a usable result. 1-4 are not.
  ok_lb <- isTRUE(res_lb$status %in% c(0L, 5L))
  note_lb <- sprintf("convergence=list(metric='marginal_kl',rule='improvement',tol=0.001) (oris_soft canonical default) requested; status=%d, iterations=%d",
                      res_lb$status, res_lb$iterations)
  row_lb <- arm_row(input_class, n, n_margins, n_categories, m_cell, max_weight,
                     "leafblower_oris_soft", as.numeric(bm_lb$median), max_error_lb,
                     max(w_lb_n), min(w_lb_n),
                     leafblower::design_effect(w_lb_n),
                     leafblower::effective_sample_size(w_lb_n),
                     res_lb$iterations, ok_lb, note_lb)
  cat(sprintf("  %-22s wall=%7.4fs status=%d max_err=%.3e max_w=%.3f n_eff=%.1f\n",
              "leafblower_oris_soft", row_lb$wall_s, res_lb$status,
              row_lb$max_error, row_lb$max_w, row_lb$n_eff))

  # --- Arm 2: survey_calibrate ---
  # Population totals derived from the SAME tgt list, scaled to n, so both
  # arms are given identical control totals (confound control). formula/
  # population as lists-of-formulas/lists-of-tables — the classical-raking
  # special case documented under ?survey::calibrate ("in the same format
  # as the input to rake"); read from the installed help, not assumed.
  formula_list <- lapply(margin_cols, function(k) stats::as.formula(paste0("~", k)))
  population_list <- lapply(margin_cols, function(k) {
    T_k <- tgt[[k]]
    data.frame(setNames(list(names(T_k)), k), Freq = as.numeric(T_k) * n)
  })
  design <- survey::svydesign(ids = ~1, weights = ~1, data = df)
  sv_call <- function() {
    survey::calibrate(design, formula_list, population_list, calfun = "raking",
                       bounds = c(0, max_weight), epsilon = 1e-3)
  }
  bm_sv  <- bench::mark(run = sv_call(), iterations = 2, check = FALSE,
                         memory = FALSE, filter_gc = FALSE)
  cal_sv <- sv_call()
  w_sv   <- stats::weights(cal_sv)
  w_sv_n <- normalize_to_n(as.numeric(w_sv), n)
  max_error_sv <- margin_max_error(w_sv_n, df, tgt)
  ok_sv <- all(is.finite(w_sv)) && max(w_sv) <= max_weight + 1e-6
  note_sv <- "epsilon=1e-3 requested (survey::calibrate, calfun='raking', bounds=c(0,max_weight))"
  # iterations is not comparable across packages for this arm.
  row_sv <- arm_row(input_class, n, n_margins, n_categories, m_cell, max_weight,
                     "survey_calibrate", as.numeric(bm_sv$median), max_error_sv,
                     max(w_sv_n), min(w_sv_n),
                     leafblower::design_effect(w_sv_n),
                     leafblower::effective_sample_size(w_sv_n),
                     NA_integer_, ok_sv, note_sv)
  cat(sprintf("  %-22s wall=%7.4fs max_err=%.3e max_w=%.3f n_eff=%.1f\n",
              "survey_calibrate", row_sv$wall_s, row_sv$max_error, row_sv$max_w, row_sv$n_eff))

  rbind(row_lb, row_sv)
}

# --- Fixture: medium_100k_5margins ---
# nj^K = 4^5 = 1024, n = 100000 -> m_cell/n far below 1: a non-degenerate,
# compression-benefiting class, deliberately NOT the K=20 uniform-random
# class RESEARCH.md Pitfall 3 warns is the wrong fixture for a positive claim.
# This shape matches the PRD's own contradictory medium-scale target
# (100K rows, 5 margins, categories per margin) so this number directly
# retires that contradiction (SC2) rather than sitting beside it.
require_single_thread_blas()

set.seed(304)
n  <- 100000L
K  <- 5L
nj <- 4L
margin_cols <- paste0("m", seq_len(K))
df <- as.data.frame(lapply(seq_len(K), function(k)
  factor(sample(letters[seq_len(nj)], n, TRUE))))
names(df) <- margin_cols

# Skewed, non-uniform per-margin target so calibration is a real adjustment,
# not a no-op. NOTE: the neighbourhood of 0.45/0.30/0.15/0.10 named in the
# plan is infeasible at max_weight=3 for K=5 INDEPENDENT margins — the
# combined multiplicative correction needed on a cell disfavoured across
# all 5 margins exceeds the bound and both arms fail to converge
# (oris_soft: status=4 budget-exhausted at max_error=3.8e-2; survey:
# "Calibration failed"). 0.40/0.28/0.18/0.14 is also too tight at n=100000
# (oris_soft hits status=5 constrained-optimum plateau at max_error=1.49e-3,
# above the 1e-3 floor task 3 gates on). 0.36/0.27/0.20/0.17 keeps the joint
# correction inside the max_weight=3 bound (leafblower still clamps weights
# at the bound on this fixture — the bound is exercised, not a no-op) while
# both arms converge (status=0) well under the 1e-3 accuracy floor. See
# 03-01-SUMMARY.md for the measured evidence at each candidate skew.
p_skew <- c(0.36, 0.27, 0.20, 0.17)
p_skew <- p_skew / sum(p_skew)
tgt <- setNames(lapply(margin_cols, function(k) setNames(p_skew, levels(df[[k]]))), margin_cols)

results <- run_input_class("medium_100k_5margins", df, tgt, max_weight = 3, n_categories = nj)

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(results, "benchmarks/results/oris_soft_vs_competitors.csv", row.names = FALSE)
cat("\nWrote benchmarks/results/oris_soft_vs_competitors.csv\n")

# --- Machine and provenance capture (SC1: the figure must name its machine) ---
cpu_model <- tryCatch({
  if (file.exists("/proc/cpuinfo")) {
    lines <- readLines("/proc/cpuinfo")
    m <- grep("^model name", lines, value = TRUE)[1]
    if (!is.na(m)) trimws(sub("^model name\\s*:\\s*", "", m)) else NA_character_
  } else {
    NA_character_  # non-Linux host: degrade to NA rather than erroring
  }
}, error = function(e) NA_character_)

si <- sessionInfo()
competitor_pkgs <- c("survey")
env_lines <- c(
  sprintf("R version: %s", R.version.string),
  sprintf("Platform: %s", R.version$platform),
  sprintf("BLAS: %s", si$BLAS),
  sprintf("LAPACK: %s", si$LAPACK),
  sprintf("CPU model: %s", if (is.na(cpu_model)) "NA (non-Linux host)" else cpu_model),
  sprintf("OMP_NUM_THREADS: %s", Sys.getenv("OMP_NUM_THREADS")),
  sprintf("OPENBLAS_NUM_THREADS: %s", Sys.getenv("OPENBLAS_NUM_THREADS")),
  sprintf("MKL_NUM_THREADS: %s", Sys.getenv("MKL_NUM_THREADS")),
  sprintf("leafblower: %s", as.character(utils::packageVersion("leafblower"))),
  vapply(competitor_pkgs, function(p)
    sprintf("%s: %s", p, as.character(utils::packageVersion(p))), character(1))
)
writeLines(env_lines, "benchmarks/results/oris_soft_vs_competitors_env.txt")
cat("Wrote benchmarks/results/oris_soft_vs_competitors_env.txt\n")
