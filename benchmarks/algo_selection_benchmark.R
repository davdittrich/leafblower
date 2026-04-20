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

# ── make_bench_data ───────────────────────────────────────────────────────────
# Generates n-row survey data with K categorical margins of cats_per_margin levels.
# Population proportions: Dirichlet(1,...,1) via normalised Exp(1) draws.
# Sample proportions: population * |Normal(1, 0.1)| noise, renormalised (~10% bias).
make_bench_data <- function(n, K, cats_per_margin) {
  stopifnot(n >= 1L, K >= 1L, cats_per_margin >= 2L)
  col_names <- paste0("m", seq_len(K))
  lvl_names <- lapply(seq_len(K), function(k) paste0("c", seq_len(cats_per_margin)))

  pop_props  <- lapply(seq_len(K), function(k) {
    x <- rexp(cats_per_margin); x / sum(x)
  })
  samp_props <- lapply(pop_props, function(p) {
    q <- p * abs(rnorm(length(p), mean = 1, sd = 0.1)); q / sum(q)
  })

  df_cols <- lapply(seq_len(K), function(k) {
    factor(sample(lvl_names[[k]], n, replace = TRUE, prob = samp_props[[k]]),
           levels = lvl_names[[k]])
  })
  df <- as.data.frame(setNames(df_cols, col_names))

  targets <- setNames(
    lapply(seq_len(K), function(k) setNames(pop_props[[k]], lvl_names[[k]])),
    col_names
  )
  list(df = df, targets = targets)
}
