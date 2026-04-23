#!/usr/bin/env Rscript
# Bayesian Level Set Estimation: ieppa (faithful) vs raking (hybrid).
# Response: y = log(t_ieppa / t_raking). Threshold: log(1.2) — ieppa wins if y <= threshold.
# Input space (3D):
#   x1 = log10(complexity) = log10(n * sum(cat_counts)) in [4, 7.7]
#   x2 = log10(tol_abs) in [-6, -3]
#   x3 = log10(prod(cat_counts) / n) in [0, 3.5]  # theoretical max compression
#
# Output: benchmarks/ieppa_vs_raking_results.rds
Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(leafblower)
  library(DiceKriging)
  library(lhs)
  source("benchmarks/plot_helpers.R")
})

.BENCH_SOURCED <- exists(".BENCH_SOURCED", envir = .GlobalEnv, inherits = FALSE)
set.seed(20260423)

X1_RANGE <- c(4.0, 7.7)
X2_RANGE <- c(-6.0, -3.0)
X3_RANGE <- c(0.0, 3.5)
THRESHOLD <- log(1.2)

# Deterministic (x1,x2,x3) → data-generation seed
bench_seed <- function(x1, x2, x3) {
  as.integer(abs(round(1e5 * sin(x1 * 13 + x2 * 17 + x3 * 19)))) %% .Machine$integer.max
}

# Generate synthetic data matching targets
generate_data <- function(x1, x2, x3) {
  seed <- bench_seed(x1, x2, x3)
  set.seed(seed)
  K <- 5
  target_product <- 10^x3 * round(10^x1 / 10)  # approximate, will be adjusted
  cats_per <- max(2L, round((10^x3 * round(10^x1 / 10))^(1/K)))
  cat_counts <- rep(cats_per, K)
  n <- round(10^x1 / sum(cat_counts))
  n <- max(n, 1000L)
  tol_abs <- 10^x2
  # Generate categorical data: uniform sampling
  cols <- lapply(seq_len(K), function(k) as.character(sample(0:(cats_per-1), n, replace = TRUE)))
  df <- as.data.frame(cols, stringsAsFactors = FALSE)
  names(df) <- paste0("v", seq_len(K))
  tgt <- lapply(seq_len(K), function(k) {
    setNames(rep(1/cats_per, cats_per), as.character(0:(cats_per-1)))
  })
  names(tgt) <- names(df)
  list(df = df, tgt = tgt, n = n, K = K, cat_counts = cat_counts, tol = tol_abs)
}

# Time a single method
time_method <- function(d, method) {
  t0 <- Sys.time()
  suppressWarnings(
    harvest(d$df, d$tgt, method = method,
            convergence = list(absolute = d$tol),
            max_iterations = 500,
            max_weight = 3)
  )
  as.numeric(Sys.time() - t0, units = "secs")
}

# Evaluate log-ratio at a single (x1,x2,x3) point
evaluate_point <- function(x1, x2, x3, runs = 3) {
  d <- generate_data(x1, x2, x3)
  t_ie <- median(replicate(runs, time_method(d, "ieppa")))
  t_rk <- median(replicate(runs, time_method(d, "raking")))
  log(t_ie / t_rk)
}

run_benchmark <- function(n_initial = 16, n_adaptive = 8) {
  # LHS initial design in [0,1]^3
  D0 <- maximinLHS(n_initial, 3)
  design <- cbind(
    X1_RANGE[1] + D0[,1] * diff(X1_RANGE),
    X2_RANGE[1] + D0[,2] * diff(X2_RANGE),
    X3_RANGE[1] + D0[,3] * diff(X3_RANGE)
  )
  colnames(design) <- c("log10_complexity", "log10_tol", "log10_compression")
  cat("Evaluating initial", n_initial, "points...\n")
  y <- numeric(n_initial)
  for (i in seq_len(n_initial)) {
    y[i] <- evaluate_point(design[i,1], design[i,2], design[i,3])
    cat(sprintf("  point %d: x=(%.2f,%.2f,%.2f) y=%.3f\n",
                i, design[i,1], design[i,2], design[i,3], y[i]))
  }
  # GP fit
  gp <- km(design = design, response = y, covtype = "matern5_2",
           nugget.estim = TRUE, control = list(trace = 0))
  # Adaptive augmentation
  for (i in seq_len(n_adaptive)) {
    # Sample candidates, pick one with highest GP predictive variance near threshold
    cand <- cbind(
      runif(1000, X1_RANGE[1], X1_RANGE[2]),
      runif(1000, X2_RANGE[1], X2_RANGE[2]),
      runif(1000, X3_RANGE[1], X3_RANGE[2])
    )
    pred <- predict(gp, cand, type = "UK", checkNames = FALSE)
    # Score: variance * indicator-near-threshold
    score <- pred$sd^2 * exp(-((pred$mean - THRESHOLD)^2) / 0.1)
    best <- cand[which.max(score), , drop = FALSE]
    y_new <- evaluate_point(best[1,1], best[1,2], best[1,3])
    design <- rbind(design, best)
    y <- c(y, y_new)
    cat(sprintf("  adaptive %d: x=(%.2f,%.2f,%.2f) y=%.3f\n",
                i, best[1,1], best[1,2], best[1,3], y_new))
    # Refit
    gp <- km(design = design, response = y, covtype = "matern5_2",
             nugget.estim = TRUE, control = list(trace = 0))
  }
  list(
    design = design, y = y, gp = gp, threshold = THRESHOLD,
    x_ranges = list(
      log10_complexity = X1_RANGE,
      log10_tol = X2_RANGE,
      log10_compression = X3_RANGE
    ),
    meta = list(seed = 20260423, n_initial = n_initial,
                n_adaptive = n_adaptive, runs_per_point = 3,
                timestamp = Sys.time())
  )
}

if (!.BENCH_SOURCED) {
  res <- run_benchmark()
  saveRDS(res, "benchmarks/ieppa_vs_raking_results.rds")
  cat("Saved benchmarks/ieppa_vs_raking_results.rds\n")
  # Plot (uses plot_helpers.R)
  # Build state object matching make_plots_3d's expected structure
  state_obj <- list(design = res$design, y = res$y, gp = res$gp,
                    bounds = list(X1_RANGE, X2_RANGE, X3_RANGE))
  make_plots_3d(state_obj, candidates = NULL, threshold = res$threshold,
                out_dir = "benchmarks", slice_dim = 3,
                slice_values = c(0, 1, 2, 3))
}
