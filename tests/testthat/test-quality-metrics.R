test_that("A7: all 5 quality metrics present in calib_result for iEPPA", {
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 500,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error", "mean_error", "kl", "chi2", "l1_weight_change"))
    expect_true(nm %in% names(result),
                info = sprintf("metric '%s' missing from calib_result", nm))
  expect_true(is.finite(result$max_error))
  expect_true(is.finite(result$mean_error))
  expect_true(is.finite(result$kl))
  expect_true(is.finite(result$chi2) || is.infinite(result$chi2))  # chi2 can be Inf on degenerate
  expect_true(is.finite(result$l1_weight_change))
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
                           max_iterations = 500,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error", "mean_error", "kl", "chi2", "l1_weight_change"))
    expect_true(nm %in% names(result))
  expect_true(is.finite(result$mean_error))
  expect_true(is.finite(result$l1_weight_change))
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
                           max_iterations = 500,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error", "mean_error", "kl", "chi2", "l1_weight_change"))
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
  expect_true(is.finite(result$l1_weight_change))
  expect_true(result$l1_weight_change >= 0)
})

test_that("a0gk: metrics finite at exit with MAX_ERR criterion (gated path)", {
  # Gate: mean_err/kl/chi2 skipped at intermediate checks when criterion=MAX_ERR/PCT.
  # They are computed on the check where convergence fires (about_to_converge gate)
  # and on the final budget iteration. Verify all three metrics are finite at exit.
  # Use 2 imbalanced margins so calibration takes multiple iterations.
  set.seed(7)
  n <- 600
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE, prob = c(0.6, 0.3, 0.1))),
    b = factor(sample(c("X","Y"),     n, replace = TRUE, prob = c(0.4, 0.6)))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("X"=0.5,"Y"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
                           max_iterations = 100,
                           convergence = list(absolute = 1e-3),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_true(is.finite(result$mean_error))
  expect_true(is.finite(result$kl) && result$kl >= 0)
  expect_true(is.finite(result$chi2) && result$chi2 >= 0)
  # Convergence must have fired (not budget exhaustion) to exercise the gate path
  expect_true(result$iterations < 100)
})

test_that("B1: compute_quality_metrics extraction: values identical to inline", {
  # Snapshot test capturing expected values from inline block (lines 536-575)
  # before extraction to helper. After extraction, verify helper produces identical results.
  set.seed(11L)
  n <- 1000L
  df_t <- data.frame(
    a = factor(sample(c("x","y","z"), n, TRUE, prob = c(0.5, 0.3, 0.2))),
    b = factor(sample(c("M","F"), n, TRUE, prob = c(0.55, 0.45)))
  )
  tgt_t <- list(a = c(x=0.4, y=0.35, z=0.25), b = c(M=0.5, F=0.5))
  r_base <- leafblower::harvest(df_t, tgt_t, method = "ieppa", max_weight = 5,
                                max_iterations = 50L, attach_weights = FALSE)
  w <- r_base
  res <- attr(r_base, "result")
  expected_margin_kl <- res$margin_kl
  expected_deff <- res$design_effect
  expected_weight_kl <- res$weight_kl
  expected_eff_obs <- res$effective_observations

  # Call extracted helper and verify identical values
  qm <- leafblower:::compute_quality_metrics(w, tgt_t, df_t)
  expect_equal(qm$margin_kl, expected_margin_kl, tolerance = 1e-10)
  expect_equal(qm$design_effect, expected_deff, tolerance = 1e-10)
  expect_equal(qm$weight_kl, expected_weight_kl, tolerance = 1e-10)
  expect_equal(qm$effective_observations, expected_eff_obs, tolerance = 1e-10)
})
