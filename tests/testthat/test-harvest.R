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
