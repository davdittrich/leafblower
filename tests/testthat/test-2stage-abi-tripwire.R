test_that("hierarchical = NULL takes early-out: hierarchical_levels_used == 0", {
  set.seed(42)
  n   <- 200L
  df  <- data.frame(
    age = sample(c("young", "old"), n, replace = TRUE, prob = c(0.6, 0.4))
  )
  tgt <- list(age = c(young = 0.5, old = 0.5))

  w <- leafblower::harvest(df, tgt, method = "raking", hierarchical = NULL)
  res <- attr(w, "result")

  # ABI tripwire: hierarchical disabled path must zero-init all diagnostic fields.
  expect_equal(res$hierarchical_levels_used, 0L)
  expect_equal(res$n_cells_total,          0L)
  expect_equal(res$n_cells_skipped,        0L)
  expect_equal(res$n_cells_inherited,      0L)
  expect_equal(res$outer_iterations_used,  0L)
  expect_equal(res$outer_residual_final,   0.0)
})
