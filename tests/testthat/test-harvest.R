test_that("method='rake' emits warning about L-BFGS-B", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a","b"), 200, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))
  expect_warning(harvest(df, tgt, method="rake"), regexp = "L-BFGS-B")
})

test_that("method='nr' emits warning about L-BFGS-B", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a","b"), 200, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))
  expect_warning(harvest(df, tgt, method="nr"), regexp = "L-BFGS-B")
})

test_that("auto-routing selects iEPPA for large complexity", {
  set.seed(1)
  n   <- 200000L
  df  <- data.frame(x = factor(sample(c("a","b","c"), n, replace=TRUE)))
  tgt <- list(x = c(a=0.33, b=0.34, c=0.33))
  # complexity = 200000 * 3 = 600000 > 500000 → iEPPA
  result <- harvest(df, tgt, method="auto")
  expect_identical(attr(result, "algorithm"), "ieppa")
})

test_that("auto-routing selects lbfgsb for small unconstrained problems", {
  set.seed(1)
  n   <- 1000L
  df  <- data.frame(x = factor(sample(c("a","b"), n, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))
  # complexity = 1000 * 2 = 2000 << 500000; max_weight=Inf (unconstrained); min_weight=0
  # New routing: L-BFGS-B only for unconstrained problems (max_weight=Inf)
  result <- harvest(df, tgt, method="auto", max_weight=Inf)
  expect_identical(attr(result, "algorithm"), "lbfgsb")
})
