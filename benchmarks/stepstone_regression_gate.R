#!/usr/bin/env Rscript
# benchmarks/stepstone_regression_gate.R
#
# Stepstone regression gate: single-stage + 2-stage paths for
# {raking, sinkhorn, greenkhorn} on stepstone_small (N~10k, K=9).
#
# Usage:
#   Rscript benchmarks/stepstone_regression_gate.R          # generate/verify
#   Rscript benchmarks/stepstone_regression_gate.R --gate   # CI mode: exit 1 on regression
#
# First-run (no fixture):  generates benchmarks/fixtures/stepstone_2stage_reference.rds
#                           and    benchmarks/fixtures/stepstone_1stage_reference.rds
# Subsequent runs:          compares against those fixtures; exits non-zero on failure.
#
# Single-stage regression guard: L∞ ≤ 1e-8 vs golden reference.
# 2-stage regression guard:
#   - outer_iterations_used  ≤ 10 (iteration budget, always satisfied by construction)
#   - outer_residual_final   ≤ fixture baseline × (1 + RESID_SLACK)  [regression]
#
# Note: the stepstone 9-margin overlapping system structurally cannot reach
# outer_tol=1e-4 in the outer loop (see benchmark comments re: cyclic IPF stall).
# The outer residual stabilises at ~0.04–0.10 on this data.  The gate records
# the golden residual on first run and checks for regression (not absolute bound).

suppressPackageStartupMessages({
  library(arrow)
  library(leafblower)
})

# ── CLI args ─────────────────────────────────────────────────────────────────
args      <- commandArgs(trailingOnly = TRUE)
gate_mode <- "--gate" %in% args

FIXTURE_1STAGE <- "benchmarks/fixtures/stepstone_1stage_reference.rds"
FIXTURE_2STAGE <- "benchmarks/fixtures/stepstone_2stage_reference.rds"

METHODS  <- c("raking", "sinkhorn", "greenkhorn")
LINF_TOL <- 1e-8     # single-stage weight regression tolerance
RESID_SLACK <- 0.05  # allow ≤ 5% relative increase in outer residual

# ── Load stepstone_small ─────────────────────────────────────────────────────
PARQUET <- "tests/testthat/fixtures/stepstone_small.parquet"
TARGETS <- "tests/testthat/fixtures/stepstone_small_targets.rds"

if (!file.exists(PARQUET))
  stop("stepstone_small.parquet not found. Run benchmarks/make_stepstone_small_fixture.R first.")
if (!file.exists(TARGETS))
  stop("stepstone_small_targets.rds not found. Run benchmarks/make_stepstone_small_fixture.R first.")

df  <- arrow::read_parquet(PARQUET)
tgt <- readRDS(TARGETS)

# Ensure factor columns for all margin variables.
for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

cat(sprintf("stepstone_small: n=%d  K=%d\n", nrow(df), length(tgt)))
cat(sprintf("margins: %s\n\n", paste(names(tgt), collapse=", ")))

# ── 2-stage config ───────────────────────────────────────────────────────────
# coarse_margins = first 3 of stepstone_small's margins (K=9).
K           <- length(tgt)
N_COARSE    <- 3L
coarse_mask <- c(rep(1L, N_COARSE), rep(0L, K - N_COARSE))

hier_cfg <- list(
  coarse_mask      = coarse_mask,
  min_cell_n       = 30L,
  mode             = 0L,      # 0L = "refine"
  outer_tol        = 1e-4,
  outer_iterations = 10L
)

cat(sprintf("2-stage coarse margins (%d): %s\n\n",
  N_COARSE, paste(names(tgt)[coarse_mask == 1L], collapse=", ")))

# ── Run single-stage (hierarchical = NULL) ────────────────────────────────────
cat("=== Single-stage (hierarchical = NULL) ===\n")
results_1stage <- setNames(vector("list", length(METHODS)), METHODS)
for (meth in METHODS) {
  t0 <- proc.time()[["elapsed"]]
  w  <- suppressWarnings(
          harvest(df, tgt, method = meth, hierarchical = NULL,
                  max_iterations = 500L, verbose = 0L, attach_weights = FALSE))
  elapsed <- proc.time()[["elapsed"]] - t0
  res     <- attr(w, "result")
  cat(sprintf("  %-12s  iters=%4d  max_err=%.3e  wall=%.2fs\n",
    meth, res$iterations, res$max_error, elapsed))
  results_1stage[[meth]] <- as.numeric(w)
}

# ── Run 2-stage (hierarchical = hier_cfg) ─────────────────────────────────────
cat("\n=== 2-stage (hierarchical=list(coarse_mask, mode=0, outer_tol=1e-4, outer_iterations=10)) ===\n")
results_2stage <- setNames(vector("list", length(METHODS)), METHODS)
for (meth in METHODS) {
  t0  <- proc.time()[["elapsed"]]
  w   <- suppressWarnings(
           harvest(df, tgt, method = meth, hierarchical = hier_cfg,
                   max_iterations = 500L, verbose = 0L, attach_weights = FALSE))
  elapsed <- proc.time()[["elapsed"]] - t0
  res     <- attr(w, "result")
  outer_iters <- res$outer_iterations_used
  outer_resid <- res$outer_residual_final
  cat(sprintf("  %-12s  outer_iters=%2d  outer_resid=%.3e  wall=%.2fs\n",
    meth, outer_iters, outer_resid, elapsed))
  results_2stage[[meth]] <- list(
    weights               = as.numeric(w),
    outer_iterations_used = outer_iters,
    outer_residual_final  = outer_resid
  )
}

# ── First-run mode: write fixtures ────────────────────────────────────────────
first_run_1stage <- !file.exists(FIXTURE_1STAGE)
first_run_2stage <- !file.exists(FIXTURE_2STAGE)

if (first_run_1stage) {
  ref1 <- list(
    git_sha      = tryCatch(system("git rev-parse --short HEAD", intern=TRUE), error=function(e) "unknown"),
    captured_at  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    K            = K,
    n            = nrow(df),
    n_coarse     = N_COARSE,
    coarse_names = names(tgt)[coarse_mask == 1L],
    weights      = results_1stage
  )
  saveRDS(ref1, FIXTURE_1STAGE)
  cat(sprintf("\n[FIRST RUN] Wrote single-stage reference: %s\n", FIXTURE_1STAGE))
}

if (first_run_2stage) {
  ref2 <- list(
    git_sha      = tryCatch(system("git rev-parse --short HEAD", intern=TRUE), error=function(e) "unknown"),
    captured_at  = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    K            = K,
    n            = nrow(df),
    n_coarse     = N_COARSE,
    coarse_names = names(tgt)[coarse_mask == 1L],
    hier_cfg     = hier_cfg,
    results      = results_2stage
  )
  saveRDS(ref2, FIXTURE_2STAGE)
  cat(sprintf("[FIRST RUN] Wrote 2-stage reference: %s\n", FIXTURE_2STAGE))
}

# ── Compare runs ──────────────────────────────────────────────────────────────
failures <- character(0)

# Single-stage: L∞ ≤ 1e-8 vs reference
cat("\n=== Regression checks ===\n")
l_inf_max <- 0.0

if (!first_run_1stage) {
  ref1 <- readRDS(FIXTURE_1STAGE)
  cat(sprintf("Single-stage reference: git_sha=%s  captured=%s\n",
    ref1$git_sha, ref1$captured_at))
  for (meth in METHODS) {
    l_inf <- max(abs(results_1stage[[meth]] - ref1$weights[[meth]]))
    l_inf_max <- max(l_inf_max, l_inf)
    pass  <- l_inf <= LINF_TOL
    cat(sprintf("  1-stage %-12s  L∞=%.3e  %s\n",
      meth, l_inf, if (pass) "PASS" else "FAIL"))
    if (!pass)
      failures <- c(failures,
        sprintf("1-stage %s: L∞=%.3e > tol=%.3e", meth, l_inf, LINF_TOL))
  }
} else {
  cat("Single-stage: first run — no comparison (reference written).\n")
}

# 2-stage: iteration budget check + regression check vs fixture
OUTER_ITER_BUDGET <- 10L
outer_iter_max    <- 0L
outer_resid_max   <- 0.0

cat("\n")
for (meth in METHODS) {
  r   <- results_2stage[[meth]]
  oi  <- r$outer_iterations_used
  or_ <- r$outer_residual_final
  outer_iter_max  <- max(outer_iter_max, oi)
  outer_resid_max <- max(outer_resid_max, or_)
  # Budget: iterations must be ≤ the requested cap (always true; verifies wiring).
  ok_iter <- oi <= OUTER_ITER_BUDGET
  cat(sprintf("  2-stage %-12s  outer_iters=%2d %s  outer_resid=%.3e\n",
    meth, oi, if (ok_iter) "" else "[FAIL-iter-budget]", or_))
  if (!ok_iter)
    failures <- c(failures,
      sprintf("2-stage %s: outer_iters=%d exceeded budget %d", meth, oi, OUTER_ITER_BUDGET))
}

# 2-stage weight regression vs fixture (rtol=1e-12) + residual non-regression
WEIGHT_RTOL <- 1e-12

if (!first_run_2stage) {
  ref2 <- readRDS(FIXTURE_2STAGE)
  cat(sprintf("\n2-stage regression vs fixture (weight L∞ tol=%.0e, resid slack=%.0f%%):\n",
    WEIGHT_RTOL, RESID_SLACK * 100))
  for (meth in METHODS) {
    w_cur   <- results_2stage[[meth]]$weights
    w_ref   <- ref2$results[[meth]]$weights
    l_inf_w <- max(abs(w_cur - w_ref))

    or_cur  <- results_2stage[[meth]]$outer_residual_final
    or_ref  <- ref2$results[[meth]]$outer_residual_final
    resid_limit <- or_ref * (1.0 + RESID_SLACK)
    ok_w    <- l_inf_w <= WEIGHT_RTOL
    ok_r    <- or_cur  <= resid_limit

    cat(sprintf("  %-12s  wt L∞=%.3e %s  resid %.3e (ref %.3e) %s\n",
      meth,
      l_inf_w, if (ok_w) "PASS" else "FAIL",
      or_cur, or_ref,   if (ok_r) "PASS" else "REGR"))
    if (!ok_w)
      failures <- c(failures,
        sprintf("2-stage weight fidelity %s: L∞=%.3e > %.3e", meth, l_inf_w, WEIGHT_RTOL))
    if (!ok_r)
      failures <- c(failures,
        sprintf("2-stage residual regression %s: %.3e > ref %.3e + %.0f%% slack",
          meth, or_cur, or_ref, RESID_SLACK * 100))
  }
} else {
  cat("2-stage regression vs fixture: first run — no comparison (reference written).\n")
}

# ── Summary ────────────────────────────────────────────────────────────────────
cat("\n=== Gate summary ===\n")
if (length(failures) == 0L) {
  cat("ALL CHECKS PASSED\n")
  cat(sprintf("  single-stage L∞_max     = %.3e (target ≤ %.0e)\n", l_inf_max, LINF_TOL))
  cat(sprintf("  2-stage outer_iter_max  = %d (budget = %d)\n", outer_iter_max, OUTER_ITER_BUDGET))
  cat(sprintf("  2-stage outer_resid_max = %.3e\n", outer_resid_max))
  if (gate_mode) quit(status = 0L)
} else {
  cat("FAILURES:\n")
  for (f in failures) cat(sprintf("  - %s\n", f))
  if (gate_mode) {
    cat("\n[GATE] Regression detected — exiting with status 1.\n")
    quit(status = 1L)
  } else {
    warning("Regression gate: ", length(failures), " failure(s) detected.")
  }
}
