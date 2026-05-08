test_that("bounds_mode='banana' (unknown) stops with error from match.arg", {
  df  <- data.frame(x = factor(c("a", "b", "a")))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  expect_error(
    harvest(df, tgt, bounds_mode = "banana", convergence = list(absolute = 1e-6)),
    regexp = "'arg' should be one of"
  )
})

test_that("bounds_mode=NA stops with error from match.arg", {
  df  <- data.frame(x = factor(c("a", "b", "a")))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  expect_error(
    harvest(df, tgt, bounds_mode = NA, convergence = list(absolute = 1e-6)),
    regexp = "'arg' must be NULL or a character vector"
  )
})

test_that("bounds_mode='cell' succeeds (default mode)", {
  set.seed(1)
  df  <- data.frame(x = factor(c("a", "b", "a")))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  suppressWarnings(
    harvest(df, tgt, bounds_mode = "cell", convergence = list(absolute = 1e-6))
  )
  expect_true(TRUE)  # placeholder: just check no error
})

test_that("bounds_mode='unit' succeeds (strict per-obs mode)", {
  set.seed(1)
  df  <- data.frame(x = factor(c("a", "b", "a")))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  suppressWarnings(
    harvest(df, tgt, bounds_mode = "unit", convergence = list(absolute = 1e-6))
  )
  expect_true(TRUE)  # placeholder: just check no error
})
