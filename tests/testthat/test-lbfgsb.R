test_that("L-BFGS-B converges on 3-margin no-bounds case", {
  set.seed(42)
  n <- 50000L
  age <- sample(c("18-34","35-54","55+"), n, replace=TRUE, prob=c(0.35,0.40,0.25))
  sex <- sample(c("M","F"), n, replace=TRUE, prob=c(0.52,0.48))
  edu <- sample(c("HS","College","Grad"), n, replace=TRUE, prob=c(0.40,0.40,0.20))
  df  <- data.frame(age=factor(age), sex=factor(sex), edu=factor(edu))
  tgt <- list(
    age = c("18-34"=0.30, "35-54"=0.45, "55+"=0.25),
    sex = c(M=0.50, F=0.50),
    edu = c(HS=0.35, College=0.45, Grad=0.20)
  )
  # RED: harvest() not yet implemented
  result <- harvest(df, tgt, method="lbfgsb")
  expect_s3_class(result, "data.frame")
  expect_true("weights" %in% names(result))
  expect_lt(abs(mean(result$weights) - 1.0), 1e-8)
})

test_that("L-BFGS-B max_weight bound respected", {
  set.seed(1)
  n <- 10000L
  x <- sample(c("a","b"), n, replace=TRUE, prob=c(0.9,0.1))
  df <- data.frame(x=factor(x))
  tgt <- list(x=c(a=0.5, b=0.5))
  result <- harvest(df, tgt, method="lbfgsb", max_weight=5)
  expect_true(max(result$weights) <= 5.0 + 1e-10)
})

test_that("near-one max_weight rejected with informative error", {
  set.seed(1)
  df <- data.frame(x = factor(sample(c("a","b"), 200, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))
  # max_weight=1.0 (exact) — denominator (U-1)=0
  expect_error(harvest(df, tgt, method="lbfgsb", max_weight=1.0),
               regexp="logit link undefined")
  # max_weight=1+5e-7 (within eps=1e-6) — catastrophic cancellation
  expect_error(harvest(df, tgt, method="lbfgsb", max_weight=1.0 + 5e-7),
               regexp="logit link undefined")
  # max_weight=1-5e-7 (within eps=1e-6) — catastrophic cancellation
  expect_error(harvest(df, tgt, method="lbfgsb", max_weight=1.0 - 5e-7),
               regexp="logit link undefined")
  # max_weight=2.0 (outside eps) — valid, should NOT error
  expect_no_error(suppressWarnings(harvest(df, tgt, method="lbfgsb", max_weight=2.0)))
})
