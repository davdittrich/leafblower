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

# ── time_cell ─────────────────────────────────────────────────────────────────
# Times iEPPA vs L-BFGS-B at one (log_complexity, log_tol) point.
# Returns log(median_t_iEPPA / median_t_LBFGSB):
#   positive → iEPPA is slower → L-BFGS-B wins.
# Threshold: log(1.2) ≈ 0.182.
#
# n derivation: cats_per_margin fixed by complexity tercile, n = round(complexity/(K*cats)).
# Both solvers run with max_weight=Inf (unconstrained, exponential link).
# seed_extra: K-specific seed offset for K-stability sweeps (default 0L = main sweep)
#   K=3 stability: pass seed_extra = 3L * 10000000L
#   K=18 stability: pass seed_extra = 18L * 10000000L
time_cell <- function(log_complexity, log_tol, K = 9L, seed_extra = 0L) {
  K <- as.integer(K)
  cats_per_margin <- if (log_complexity <= 5.5) 4L else if (log_complexity <= 6.5) 8L else 16L
  n <- max(50L, 2L * K * cats_per_margin,
         as.integer(round(10^log_complexity / (K * cats_per_margin))))

  # Round-trip sanity check
  actual_log_c <- log10(n * K * cats_per_margin)
  if (abs(actual_log_c - log_complexity) > 0.2)
    warning(sprintf("time_cell: complexity round-trip %.2f log-units (lc=%.2f K=%d cats=%d n=%d)",
                    abs(actual_log_c - log_complexity), log_complexity, K, cats_per_margin, n))

  # Generate data; retry with progressively larger n if any margin has empty cells.
  # Empty cells make harvest() infeasible regardless of solver.
  seed_base <- bench_seed(log_complexity, log_tol) + seed_extra
  bd <- NULL
  n_try <- n
  for (.attempt in seq_len(5L)) {
    set.seed(seed_base + .attempt - 1L)
    cand <- make_bench_data(n_try, K, cats_per_margin)
    all_filled <- all(sapply(cand$df, function(col) min(table(col)) > 0L))
    if (all_filled) { bd <- cand; break }
    n_try <- as.integer(ceiling(n_try * 1.5))
  }
  if (is.null(bd)) stop("time_cell: could not generate feasible data after 5 attempts")
  conv <- list(absolute = 10^log_tol)

  time_algo <- function(method) {
    # 2 warmup runs (discarded)
    for (i in seq_len(2L))
      suppressWarnings(invisible(leafblower::harvest(
        bd$df, bd$targets, method = method,
        max_weight = Inf, min_weight = 0, convergence = conv, max_iterations = 500L)))
    # 5 timed runs
    median(replicate(5L, {
      t0 <- proc.time()[["elapsed"]]
      suppressWarnings(invisible(leafblower::harvest(
        bd$df, bd$targets, method = method,
        max_weight = Inf, min_weight = 0, convergence = conv, max_iterations = 500L)))
      proc.time()[["elapsed"]] - t0
    }))
  }

  t_ieppa  <- time_algo("ieppa")
  t_lbfgsb <- time_algo("lbfgsb")
  # Floor at 0.1ms to guard against proc.time() resolution yielding exact zeros.
  t_min <- 1e-4
  log(max(t_ieppa, t_min) / max(t_lbfgsb, t_min))
}

# ── fit_gp ────────────────────────────────────────────────────────────────────
# Fits a Matérn-5/2 GP to timing observations.
# design_mat: n_obs × 2 matrix (col 1 = log_complexity, col 2 = log_tol).
# y: numeric vector of log(t_iEPPA / t_LBFGSB).
fit_gp <- function(design_mat, y) {
  stopifnot(is.matrix(design_mat), nrow(design_mat) == length(y), ncol(design_mat) == 2L)
  DiceKriging::km(
    formula      = ~1,
    design       = as.data.frame(design_mat),
    response     = y,
    covtype      = "matern5_2",
    nugget.estim = TRUE,
    nugget       = 1e-4,          # lower bound: prevents degenerate fit on small n
    control      = list(trace = FALSE)
  )
}

# ── straddle_next ─────────────────────────────────────────────────────────────
# Straddle acquisition (Bryan et al. 2005): picks the candidate maximising
#   a(x) = -|μ(x) − threshold| + κ·σ(x)
# Pulls samples toward the contour (low |μ − threshold|) and uncertain regions (high σ).
straddle_next <- function(gp_model, candidates, threshold, kappa = 2) {
  pred <- DiceKriging::predict(gp_model,
                               newdata    = as.data.frame(candidates),
                               type       = "UK",
                               checkNames = FALSE)
  a   <- -abs(pred$mean - threshold) + kappa * pred$sd
  candidates[which.max(a), , drop = FALSE]
}

# ── classified_fraction ───────────────────────────────────────────────────────
# Fraction of candidates classified with ≥conf confidence as above or below threshold.
# Termination fires when this reaches 0.90.
classified_fraction <- function(gp_model, candidates, threshold, conf = 0.95) {
  pred    <- DiceKriging::predict(gp_model,
                                  newdata    = as.data.frame(candidates),
                                  type       = "UK",
                                  checkNames = FALSE)
  p_above <- pnorm(threshold, mean = pred$mean, sd = pred$sd, lower.tail = FALSE)
  mean(p_above > conf | p_above < (1 - conf))
}

# ── save_checkpoint ───────────────────────────────────────────────────────────
# Atomic checkpoint: write to .tmp, then rename to final path.
# state fields: design (n×2 matrix), y (numeric), gp (km or NULL),
#               iter (integer), classified (numeric), bounds (list).
save_checkpoint <- function(state, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(state, tmp)
  file.rename(tmp, path)   # atomic on same filesystem (benchmarks/ → benchmarks/)
  invisible(path)
}

# ── load_checkpoint ───────────────────────────────────────────────────────────
# Returns NULL if no checkpoint found; otherwise returns the saved state.
load_checkpoint <- function(path) {
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}
