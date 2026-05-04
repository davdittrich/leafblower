# tests/testthat/test-safety.R
library(leafblower)

.mini_harvest <- function(n = 20, method = "raking") {
  set.seed(1)
  df <- data.frame(
    sex  = sample(c("M", "F"), n, replace = TRUE),
    wt   = rep(1.0, n)
  )
  target <- list(sex = c(M = 0.5, F = 0.5))
  harvest(df, target = target, weight_column = "wt", method = method)
}

test_that("B4: solver dispatch returns without crashing for each method", {
  for (m in c("raking", "ieppa", "auto")) {
    expect_no_error(.mini_harvest(method = m))
  }
})

test_that("B7: start_weights length mismatch produces error", {
  set.seed(2)
  n  <- 30L
  df <- data.frame(sex = sample(c("M", "F"), n, replace = TRUE))
  tgt <- list(sex = c(M = 0.5, F = 0.5))
  df[["sex"]] <- factor(df[["sex"]], levels = c("M","F"))
  bad_sw <- rep(1.0, n - 5L)  # wrong length: 25 instead of 30
  expect_error(
    harvest(df, tgt, start_weights = bad_sw),
    regexp = "start_weights"
  )
})
