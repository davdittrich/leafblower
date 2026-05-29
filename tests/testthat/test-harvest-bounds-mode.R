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
  n <- 200
  df  <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  r <- suppressWarnings(
    harvest(df, tgt, bounds_mode = "cell", convergence = list(absolute = 1e-6))
  )
  expect_equal(attr(r, "result")$status, 0L, info = "cell mode must converge")
  expect_lt(abs(sum(r$weights) - nrow(df)), 1e-6,
            label = "sum(weights) == n")
})

test_that("bounds_mode='unit' succeeds (strict per-obs mode)", {
  set.seed(1)
  n <- 200
  max_weight <- 3.0
  min_weight <- 0.2
  df  <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  r <- suppressWarnings(
    harvest(df, tgt, bounds_mode = "unit",
            max_weight = max_weight, min_weight = min_weight,
            convergence = list(absolute = 1e-6))
  )
  expect_equal(attr(r, "result")$status, 0L, info = "unit mode must converge")
  expect_lte(max(r$weights), max_weight + 1e-9, label = "max(weights) <= max_weight")
  expect_gte(min(r$weights), min_weight - 1e-9, label = "min(weights) >= min_weight")
  expect_lt(abs(sum(r$weights) - nrow(df)), 1e-6, label = "sum(weights) == n")
})
