test_that("rk_calibrate honours caller-supplied max_weight (not overwritten by defaults)", {
  # This regression test verifies that rk_params_init is ONLY called for
  # default initialization when params == NULL. If someone accidentally changes
  # the C code to call rk_params_init(params), this test will fail because
  # max_weight would be reset to the default (5.0) instead of honouring 2.0.

  set.seed(42)
  n   <- 100L
  df  <- data.frame(x = factor(rep("a", n)))
  tgt <- list(x = c(a = 1.0))

  # Pass max_weight=2.0. If rk_params_init were called on the supplied struct,
  # it would reset to the default (5.0) and weights could exceed 2.0.
  result <- harvest(df, tgt, method = "raking", max_weight = 2.0,
                    convergence = list(absolute = 1e-6))
  weights <- result$weights

  # Check that max_weight=2.0 was honoured: all weights <= 2.0 (with tolerance)
  expect_lte(max(weights), 2.0 + 1e-9)
})

test_that("rk_calibrate honours caller-supplied min_weight", {
  # Companion test: verify min_weight is also honoured when caller-supplied.
  set.seed(42)
  n   <- 100L
  df  <- data.frame(x = factor(rep("a", n)))
  tgt <- list(x = c(a = 1.0))

  # Pass min_weight=0.5. If rk_params_init were called on the supplied struct,
  # it would reset to the default (0.0) and weights could drop below 0.5.
  result <- harvest(df, tgt, method = "raking", min_weight = 0.5,
                    convergence = list(absolute = 1e-6))
  weights <- result$weights

  # Check that min_weight=0.5 was honoured: all weights >= 0.5 (with tolerance)
  expect_gte(min(weights), 0.5 - 1e-9)
})
