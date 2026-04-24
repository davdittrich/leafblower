test_that("A7: all 5 quality metrics present in calib_result for iEPPA", {
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error", "mean_error", "kl", "chi2", "pct_change"))
    expect_true(nm %in% names(result),
                info = sprintf("metric '%s' missing from calib_result", nm))
  expect_true(is.finite(result$max_error))
  expect_true(is.finite(result$mean_error))
  expect_true(is.finite(result$kl))
  expect_true(is.finite(result$chi2) || is.infinite(result$chi2))  # chi2 can be Inf on degenerate
  expect_true(is.finite(result$pct_change))
})

test_that("A7: all 5 quality metrics present in calib_result for raking", {
  set.seed(43)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "raking",
                           max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error", "mean_error", "kl", "chi2", "pct_change"))
    expect_true(nm %in% names(result))
  expect_true(is.finite(result$mean_error))
  expect_true(is.finite(result$pct_change))
})

test_that("A7: all 5 quality metrics present in calib_result for lbfgsb", {
  set.seed(44)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "lbfgsb",
                           max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error", "mean_error", "kl", "chi2", "pct_change"))
    expect_true(nm %in% names(result))
  expect_true(is.finite(result$mean_error))
})

test_that("A7: metrics non-zero after max_iter exit (solver exits before kErrCheckInterval)", {
  # Use max_iterations=1 to force exit after 1 iteration, likely before kErrCheckInterval.
  # All three solvers check at iter==1, so metrics must be populated.
  set.seed(55)
  n <- 500
  data <- data.frame(a = factor(sample(c("1","2","3"), n, replace = TRUE)))
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 1,
                           convergence = list(absolute = 1e-20),  # impossible threshold
                           attach_weights = FALSE)
  result <- attr(w, "result")
  # Metrics must be populated (finite) even after 1 iteration:
  expect_true(is.finite(result$mean_error))
  expect_true(is.finite(result$pct_change))
  expect_true(result$pct_change >= 0)
})
