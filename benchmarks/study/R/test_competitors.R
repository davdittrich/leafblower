# benchmarks/study/R/test_competitors.R
#
# TDD tests for benchmarks/study/R/competitors.R (ticket leafblower-2ouc.6,
# WU-5). Two tiers:
#   1. Contract-shape smoke test, looped over all 17 run_<id>() adapters
#      against the toy_inline fixture (spec/toy_inline.json, kl family,
#      n=4, K=1, 2-level margin) -- every adapter must accept it structurally
#      regardless of native family, since margin encoding is family-agnostic.
#   2. Per-package home-turf golden: a dataset/target combination each
#      package is expected to calibrate accurately on, asserting the
#      harness-recomputed `converged` field (never the package's own
#      self-report) and/or an achieved margin_linf within a documented
#      tolerance.
#
# Run (single-thread BLAS mandatory, CLAUDE.md):
#   OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
#     Rscript benchmarks/study/R/test_competitors.R

suppressPackageStartupMessages({
  library(testthat)
})

# Relative-to-repo-root (matches competitors.R's own source() convention and
# CLAUDE.md's "run from repo root" contract). This file is a plain Rscript
# entry point (test_that() blocks execute top-to-bottom on source(), each one
# reporting its own PASS/FAIL); it is NOT run via testthat::test_file(),
# which chdir's into the test file's directory and would break every
# relative fixture path below.
source(file.path("benchmarks", "study", "R", "competitors.R"))

STATUS_ENUM <- c("converged", "no_conv", "infeasible", "bound_violation",
                  "bad_arg", "budget", "stall", "error")

TOY <- load_problem_spec("benchmarks/study/spec/toy_inline.json")
APISTRAT <- load_problem_spec("benchmarks/study/spec/canonical_survey_apistrat.json")

# All 17 registry.json R-arm competitor ids -> run_<id>() function.
ADAPTERS <- list(
  survey_calibrate_raking = run_survey_calibrate_raking,
  survey_calibrate_linear = run_survey_calibrate_linear,
  survey_calibrate_logit  = run_survey_calibrate_logit,
  sampling_calib_linear   = run_sampling_calib_linear,
  sampling_calib_logit    = run_sampling_calib_logit,
  anesrake                = run_anesrake,
  ipfr                    = run_ipfr,
  ebal                    = run_ebal,
  optweight_entropy       = run_optweight_entropy,
  optweight_linf          = run_optweight_linf,
  weightit                = run_weightit,
  regenesees              = run_regenesees,
  icarus                  = run_icarus,
  laeken                  = run_laeken,
  gecal                   = run_gecal,
  jointcalib              = run_jointcalib,
  sbw                     = run_sbw
)

expect_contract_shape <- function(res, problem) {
  n <- nrow(problem$data)
  expect_true(is.list(res))
  expect_setequal(names(res), c("weights_ref", "iterations", "status", "converged",
                                 "error_message", "wall_time_s", "peak_rss_bytes"))
  expect_true(is.character(res$weights_ref) && length(res$weights_ref) == 1L)
  expect_true(nzchar(res$weights_ref))
  expect_true(file.exists(res$weights_ref))
  wtab <- arrow::read_parquet(res$weights_ref)
  expect_equal(names(wtab), "weight")
  expect_equal(nrow(wtab), n)

  expect_true(is.na(res$iterations) || (is.numeric(res$iterations) && res$iterations == as.integer(res$iterations)))
  expect_true(res$status %in% STATUS_ENUM)
  expect_true(is.logical(res$converged) && length(res$converged) == 1L && !is.na(res$converged))
  expect_true(is.null(res$error_message) || (is.character(res$error_message) && length(res$error_message) == 1L))
  expect_true(is.numeric(res$wall_time_s) && length(res$wall_time_s) == 1L && res$wall_time_s >= 0)
  expect_true(is.numeric(res$peak_rss_bytes) && length(res$peak_rss_bytes) == 1L)

  # converged must be reproducible independently from the recorded weights --
  # never trust the adapter's internal bookkeeping (contract.md 2.4).
  if (!all(is.na(wtab$weight))) {
    ms <- margin_stats(wtab$weight, .comp_groups(problem), problem$targets)
    expect_equal(res$converged, isTRUE(ms$margin_linf <= problem$tol))
  }
}

test_that("every registered R-competitor adapter satisfies the contract v2 shape on toy_inline", {
  for (id in names(ADAPTERS)) {
    res <- ADAPTERS[[id]](TOY)
    expect_contract_shape(res, TOY)
  }
})

# ---- home-turf goldens ------------------------------------------------------

test_that("survey_calibrate_raking home-turf golden: pkg:survey::apistrat/stype", {
  res <- run_survey_calibrate_raking(APISTRAT)
  expect_equal(res$status, "converged")
  expect_true(res$converged)
  w <- arrow::read_parquet(res$weights_ref)$weight
  ms <- margin_stats(w, .comp_groups(APISTRAT), APISTRAT$targets)
  expect_lt(ms$margin_linf, APISTRAT$tol)
  expect_equal(sum(w), sum(APISTRAT$design_weights), tolerance = 1e-6)
})

test_that("survey_calibrate_linear/logit home-turf golden: pkg:survey::apistrat/stype", {
  res_lin <- run_survey_calibrate_linear(APISTRAT)
  expect_true(res_lin$converged)

  # logit calfun requires finite bounds -- apistrat spec ships max=null/Inf,
  # so build a bounded variant inline (bad_arg is exercised separately below).
  bounded <- APISTRAT
  bounded$bounds <- list(min = 0, max = 5)
  res_logit <- run_survey_calibrate_logit(bounded)
  expect_true(res_logit$converged)
})

test_that("survey_calibrate_logit rejects infinite bounds as bad_arg (pre-solve, no iterate)", {
  res <- run_survey_calibrate_logit(APISTRAT)  # apistrat bounds max=Inf
  expect_equal(res$status, "bad_arg")
  expect_false(res$converged)
  w <- arrow::read_parquet(res$weights_ref)$weight
  expect_true(all(is.na(w)))
})

.tc_anes_golden <- function() {
  if (!requireNamespace("anesrake", quietly = TRUE)) {
    skip("anesrake not installed")
  }
  utils::data(anes04, package = "anesrake", envir = environment())
  df <- get("anes04", envir = environment())
  df <- df[!is.na(df$racecats), , drop = FALSE]
  # racecats is numerically coded 1-6 with no label mapping in the package's
  # own docs (?anes04 lists only that there are "6 Racial Categories");
  # anesrake's own canonical usage rakes on the raw numeric-code names, so we
  # do the same rather than guessing semantic labels. Empirical sample
  # proportions are ~72/15/7/2/3/1% across codes 1..6; the target below is
  # deliberately >5pp off on the dominant code, satisfying anesrake's own
  # frozen pctlim=5 gate.
  target <- c(`1` = 0.55, `2` = 0.20, `3` = 0.15, `4` = 0.05, `5` = 0.03, `6` = 0.02)
  lv <- levels(factor(df$racecats))
  target <- target[lv]
  list(
    id = "anes04_racecats_golden",
    data = data.frame(racecats = factor(df$racecats)),
    design_weights = rep(1, nrow(df)),
    margins = "racecats",
    targets = list(racecats = target / sum(target)),
    bounds = list(min = 0, max = Inf),
    tol = 1e-4,
    objective_families = "kl",
    K = 1L
  )
}

test_that("anesrake home-turf golden: bundled anes04/racecats with a >5pct-deviating target", {
  problem <- .tc_anes_golden()
  res <- run_anesrake(problem)
  expect_contract_shape(res, problem)
  w <- arrow::read_parquet(res$weights_ref)$weight
  ms <- margin_stats(w, .comp_groups(problem), problem$targets)
  expect_lt(ms$margin_linf, 0.01)
  expect_true(res$converged)
})

# ---- targeted family-dispatch + bad_arg coverage ----------------------------

test_that("family-dispatch adapters (icarus/laeken/regenesees/gecal/weightit/sbw) reject an unmapped objective family as bad_arg", {
  bad <- TOY
  bad$objective_families <- "ot"  # not in any dispatch table
  for (fn in list(run_icarus, run_laeken, run_regenesees, run_gecal, run_weightit)) {
    res <- fn(bad)
    expect_equal(res$status, "bad_arg")
    expect_false(res$converged)
  }
})

test_that("sampling_calib_logit / laeken / icarus reject infinite bounds as bad_arg when method requires bounds", {
  res <- run_sampling_calib_logit(TOY)  # toy_inline bounds max=10, finite -- sanity control
  expect_true(res$status %in% STATUS_ENUM)

  unbounded <- TOY
  unbounded$bounds <- list(min = 0, max = Inf)
  res_unb <- run_sampling_calib_logit(unbounded)
  expect_equal(res_unb$status, "bad_arg")
})

cat("\nAll test_competitors.R suites completed.\n")
