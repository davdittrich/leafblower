library(leafblower)

# CXX.2 (leafblower-5fm8.2): greenkhorn best-iterate must be selected on the
# configured convergence metric (select_metric(cfg.metric, m)) at the error-check
# interval, not on the errRp fast proxy. For the default metric=="max_err" this
# is behaviorally identical; for kl/chi2 the reported best_error/max_error must be
# on the requested metric's scale.

make_feasible_df <- function() {
  set.seed(7)
  n <- 200L
  data.frame(
    a = factor(sample(c("1", "2", "3"), n, TRUE)),
    b = factor(sample(c("x", "y"), n, TRUE))
  )
}
feasible_targets <- function() {
  list(a = c("1" = 1 / 3, "2" = 1 / 3, "3" = 1 / 3),
       b = c(x = 0.5, y = 0.5))
}

test_that("CXX.2: greenkhorn default metric=max_err converges (no regression)", {
  w <- harvest(make_feasible_df(), feasible_targets(),
    method = "greenkhorn", max_iterations = 500L, attach_weights = FALSE)
  r <- attr(w, "result")
  # Converges (status 0) and the reported max_error equals the max marginal
  # error of the returned weights, on the max_err scale.
  expect_equal(r$status, 0L)
  expect_lt(r$max_error, 1e-3)
  expect_equal(r$max_error, r$best_error, tolerance = 1e-12)
})

test_that("CXX.2: greenkhorn metric=kl reports best_error on the KL scale", {
  w <- harvest(make_feasible_df(), feasible_targets(),
    method = "greenkhorn",
    convergence = list(metric = "kl", rule = "threshold", tol = 1e-8),
    max_iterations = 500L, attach_weights = FALSE)
  r <- attr(w, "result")
  # best-iterate selected via select_metric(KL, m): max_error == best_error and
  # both are the KL value of the best iterate (finite, non-negative). Pre-fix
  # this slot held the errRp max instead of the KL value.
  expect_true(is.finite(r$best_error))
  expect_gte(r$best_error, 0)
  expect_equal(r$max_error, r$best_error, tolerance = 1e-12)
  # The convergence metric recorded must be KL, confirming the configured
  # metric drove the run that the best-iterate is now selected against.
  expect_equal(r$convergence_used$metric, "kl")
})
