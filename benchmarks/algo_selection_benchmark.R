# benchmarks/algo_selection_benchmark.R
# Bayesian Level Set Estimation benchmark for iEPPA vs L-BFGS-B algorithm selection.
# See: docs/superpowers/specs/2026-04-20-algo-selection-design.md

# Enforce single-threaded execution before loading any library that might
# initialise an OpenMP thread pool. OMP reads OMP_NUM_THREADS at pool creation
# time (lazy, first parallel region). Setting it here — before library() calls —
# guarantees the env var is visible at that point, eliminating CPU-contention
# confound from the timing comparisons.
Sys.setenv(OMP_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(leafblower)
  library(DiceKriging)
  library(lhs)
  library(ggplot2)
})

# Guard: set .BENCH_SOURCED <- TRUE before source()ing this file in tests
# to prevent run_benchmark() from executing.
.BENCH_SOURCED <- exists(".BENCH_SOURCED", envir = .GlobalEnv, inherits = FALSE)

# Input space bounds
BENCH_X1_RANGE <- c(4.0, 7.7)   # log10(complexity): 10K to 50M
BENCH_X2_RANGE <- c(-6.0, -3.0) # log10(tol_abs):    1e-6 to 1e-3
BENCH_THRESHOLD <- log(1.2)      # log(1.2) ≈ 0.182; L-BFGS-B wins above this

# ── bench_seed ────────────────────────────────────────────────────────────────
# Deterministic integer seed from (log_complexity, log_tol).
# Valid for log_complexity in [4, 7.7] and log_tol in [-6, -3].
# Formula produces seeds in [4e8, 7.7e8+5999], all within 32-bit integer range.
bench_seed <- function(log_complexity, log_tol) {
  a <- as.integer(round(log_complexity * 1e4))
  b <- as.integer(round(-log_tol * 1e4)) %% 10000L
  (a * 10000L) + b
}
