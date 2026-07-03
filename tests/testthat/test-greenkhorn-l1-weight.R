# ──────────────────────────────────────────────────────────────────────────────
# CR-A7 (mxcl.7): greenkhorn + metric="l1_weight" must not spuriously converge on
# a phantom zero.
#
# CellMetrics.l1 (obs-level Σ_i|w_out−w_in|/n) is written only by raking.cpp and
# sinkhorn.cpp before their convergence checks; compute_cell_metrics does NOT
# compute it. greenkhorn drives convergence through compute_cell_metrics +
# select_metric, so metric="l1_weight" fed a phantom m.l1 == 0.0 into
# check_convergence, which instant-fired at the first check interval (iter 10,
# achieved max_err ~1.9e-4) instead of converging like metric="max_err".
#
# greenkhorn's greedy coordinate-wise structure has no per-check prev-weight
# snapshot, so — mirroring chebyshev.cpp — the fix remaps L1_WEIGHT → MAX_ERR at
# solver entry, tracking the actual objective. (chebyshev already does this;
# raking/sinkhorn populate a real l1 and are unaffected.)
# ──────────────────────────────────────────────────────────────────────────────

test_that("greenkhorn metric=l1_weight converges, not phantom-zero spurious OK (CR-A7)", {
  set.seed(707); n <- 2000L
  a <- factor(sample(c("A", "B", "C"), n, TRUE, prob = c(.6, .25, .15)))
  b <- factor(sample(c("X", "Y"),      n, TRUE, prob = c(.7, .3)))
  data   <- data.frame(a = a, b = b)
  target <- list(a = c(A = .33, B = .34, C = .33), b = c(X = .5, Y = .5))

  achieved <- function(w) {
    Wt <- sum(w)
    max(vapply(names(target), function(v)
      max(abs(tapply(w, data[[v]], sum)[names(target[[v]])] / Wt - target[[v]])),
      numeric(1)))
  }

  w_l1 <- suppressWarnings(harvest(data, target, method = "greenkhorn",
                                   convergence = list(metric = "l1_weight"),
                                   max_iterations = 500L, attach_weights = FALSE))
  w_max <- suppressWarnings(harvest(data, target, method = "greenkhorn",
                                    convergence = list(metric = "max_err"),
                                    max_iterations = 500L, attach_weights = FALSE))

  # Must NOT phantom-stop at iter 10 (~1.9e-4); must reach the real objective.
  expect_lt(achieved(w_l1), 1e-6, label = "greenkhorn l1_weight achieved")
  # Remapped to MAX_ERR, so l1_weight and max_err reach the same fixed point.
  expect_equal(achieved(w_l1), achieved(w_max), tolerance = 1e-9,
               label = "greenkhorn l1_weight == max_err (remapped)")
})
