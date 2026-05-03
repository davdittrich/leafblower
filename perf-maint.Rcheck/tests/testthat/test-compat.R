test_that("anesrake wraps harvest with rake->lbfgsb warning", {
  set.seed(42)
  n <- 100L
  df <- data.frame(x = factor(sample(c("a","b"), n, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))
  expect_warning(
    result <- anesrake(df, tgt, choosemethod="rake"),
    regexp = "not implemented"
  )
  expect_true(is.numeric(result$weights))
  expect_equal(length(result$weights), n)
})

test_that("get_current_miss returns max calibration error", {
  set.seed(1)
  n <- 200L
  df <- data.frame(x = factor(sample(c("a","b"), n, replace=TRUE, prob=c(0.7,0.3))))
  tgt <- list(x = c(a=0.5, b=0.5))
  result <- harvest(df, tgt, convergence = list(absolute = 1e-6))
  miss <- get_current_miss(df, tgt, result$weights)
  expect_true(is.numeric(miss))
  expect_true(miss >= 0)
  expect_true(miss < 1e-3)
})

test_that("weighted_pct computes weighted proportions", {
  w <- c(1, 2, 3, 4)
  x <- factor(c("a", "b", "a", "b"))
  pct <- weighted_pct(x, w)
  expect_true(is.numeric(pct))
  expect_equal(sum(pct), 1.0, tolerance = 1e-12)
  expect_equal(pct[["a"]], 0.4, tolerance = 1e-10)
  expect_equal(pct[["b"]], 0.6, tolerance = 1e-10)
})
