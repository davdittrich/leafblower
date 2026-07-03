# CR-D14 (leafblower-j7x8.14): design_effect uses 64-bit (size_t) index arithmetic
# and an n*K/n*p overflow guard. Regression lock: the index cast must not corrupt
# normal indexing (the overflow itself needs n~34M rows, verified by inspection).

test_that("design_effect indexing correct after 64-bit conversion (CR-D14)", {
  set.seed(6); n <- 500
  df <- data.frame(a = factor(sample(c("x", "y", "z"), n, TRUE)),
                   b = factor(sample(c("p", "q"), n, TRUE)))
  tg <- list(a = c(x = 1/3, y = 1/3, z = 1/3), b = c(p = .5, q = .5))
  r <- attr(suppressWarnings(harvest(df, tg, method = "raking", max_iter = 200)), "result")
  expect_true(is.finite(r$design_effect))
  expect_gt(r$design_effect, 0)
})
