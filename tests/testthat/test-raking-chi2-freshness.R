# Regression guard for leafblower-24f7:
# chi2 in the raking result must reflect the convergence-firing iterate's weights,
# not a stale value from a prior full-metrics evaluation.
#
# Proof strategy: independently recompute chi2 from the returned final weights
# using the exact Pearson formula from compute_cell_metrics (calib_dispatch.hpp:263):
#   chi2 = sum_{k,j: T_kj*n > kMetricEps} (obs_kj - T_kj*n)^2 / (T_kj*n)
# where obs_kj = sum of weights in margin k, category j, and n = nrow(df).
#
# Three scenarios tested:
#   1. Flat loop, metric="max_err", convergence at iter=1 (first check iteration).
#   2. Flat loop, metric="max_err", convergence at iter=kErrCheckInterval=10
#      (forces >=10 iters, confirms freshness on later check iter too).
#   3. SRAA path, metric="max_err", early convergence.

# Helper: recompute Pearson chi2 from weights and targets exactly as C does.
# W = n (raking normalises sum(weights) = n, and compute_cell_metrics uses st.n).
# kMetricEps = 1e-10 (calib_dispatch.hpp:243).
recompute_chi2 <- function(weights, targets, df) {
  n <- length(weights)
  kMetricEps <- 1e-10
  chi2 <- 0.0
  for (nm in names(targets)) {
    tgt_k <- targets[[nm]]
    cats   <- names(tgt_k)
    fac    <- df[[nm]]
    for (cat in cats) {
      T_j  <- tgt_k[[cat]]
      pop  <- T_j * n                          # T_kj * W, W = n
      if (pop > kMetricEps) {
        obs  <- sum(weights[fac == cat])       # obs_kj
        chi2 <- chi2 + (obs - pop)^2 / pop
      }
    }
  }
  chi2
}

# --- Scenario 1: flat loop, converges at iter=1 (metric="max_err") ---
test_that("24f7: chi2 is fresh when flat-loop MAX_ERR converges at iter=1", {
  # Perfectly balanced problem: sample exactly matches targets so iter=1 converges.
  set.seed(1L)
  n   <- 200L
  # Construct df so empirical proportions exactly equal targets — convergence in 1 iter.
  df  <- data.frame(
    x = factor(c(rep("A", 100L), rep("B", 100L)))
  )
  tgt <- list(x = c(A = 0.5, B = 0.5))

  res <- leafblower::harvest(df, tgt, method = "raking",
                             convergence = list(metric = "max_err", absolute = 1e-6),
                             accelerate = FALSE,
                             max_iterations = 500L)
  raw   <- attr(res, "result")
  iters <- raw$iterations
  # Should converge very early (iter 1 or kErrCheckInterval=10).
  expect_lte(iters, 10L)

  chi2_reported    <- raw$chi2
  chi2_independent <- recompute_chi2(res$weights, tgt, df)

  expect_equal(chi2_reported, chi2_independent, tolerance = 1e-9,
               label = sprintf("chi2 fresh at iter=%d: reported=%.6e independent=%.6e",
                               iters, chi2_reported, chi2_independent))
})

# --- Scenario 2: flat loop, forced to converge at a later check-interval iter ---
test_that("24f7: chi2 is fresh when flat-loop MAX_ERR converges at iter>1", {
  # Use a problem that requires ~10-20 iters to converge: slight imbalance.
  set.seed(2L)
  n   <- 1000L
  df  <- data.frame(
    x = factor(sample(c("A","B"), n, replace = TRUE, prob = c(0.48, 0.52))),
    y = factor(sample(c("P","Q"), n, replace = TRUE, prob = c(0.51, 0.49)))
  )
  tgt <- list(x = c(A = 0.5, B = 0.5), y = c(P = 0.5, Q = 0.5))

  res <- leafblower::harvest(df, tgt, method = "raking",
                             convergence = list(metric = "max_err", absolute = 1e-4),
                             accelerate = FALSE,
                             max_iterations = 500L)
  raw   <- attr(res, "result")
  iters <- raw$iterations

  chi2_reported    <- raw$chi2
  chi2_independent <- recompute_chi2(res$weights, tgt, df)

  expect_equal(chi2_reported, chi2_independent, tolerance = 1e-9,
               label = sprintf("chi2 fresh at iter=%d: reported=%.6e independent=%.6e",
                               iters, chi2_reported, chi2_independent))
})

# --- Scenario 3: SRAA path, early convergence ---
test_that("24f7: chi2 is fresh on SRAA path with early convergence", {
  set.seed(3L)
  n   <- 500L
  df  <- data.frame(
    x = factor(sample(c("A","B"), n, replace = TRUE, prob = c(0.48, 0.52))),
    y = factor(sample(c("P","Q"), n, replace = TRUE, prob = c(0.53, 0.47)))
  )
  tgt <- list(x = c(A = 0.5, B = 0.5), y = c(P = 0.5, Q = 0.5))

  res <- leafblower::harvest(df, tgt, method = "raking",
                             convergence = list(metric = "max_err", absolute = 1e-4),
                             accelerate = TRUE,
                             max_iterations = 500L)
  raw   <- attr(res, "result")
  iters <- raw$iterations

  chi2_reported    <- raw$chi2
  chi2_independent <- recompute_chi2(res$weights, tgt, df)

  expect_equal(chi2_reported, chi2_independent, tolerance = 1e-9,
               label = sprintf("chi2 fresh (SRAA) at iter=%d: reported=%.6e independent=%.6e",
                               iters, chi2_reported, chi2_independent))
})

# --- Scenario 4: flat loop, metric="chi2" (chi2 is the convergence criterion) ---
test_that("24f7: chi2 is fresh when flat-loop convergence metric is chi2", {
  set.seed(4L)
  n   <- 800L
  df  <- data.frame(
    x = factor(sample(c("A","B","C"), n, replace = TRUE, prob = c(0.4, 0.35, 0.25)))
  )
  tgt <- list(x = c(A = 0.33, B = 0.34, C = 0.33))

  res <- leafblower::harvest(df, tgt, method = "raking",
                             convergence = list(metric = "chi2", absolute = 1e-3),
                             accelerate = FALSE,
                             max_iterations = 500L)
  raw   <- attr(res, "result")
  iters <- raw$iterations

  chi2_reported    <- raw$chi2
  chi2_independent <- recompute_chi2(res$weights, tgt, df)

  expect_equal(chi2_reported, chi2_independent, tolerance = 1e-9,
               label = sprintf("chi2 fresh (chi2-metric) at iter=%d: reported=%.6e independent=%.6e",
                               iters, chi2_reported, chi2_independent))
})
