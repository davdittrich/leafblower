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
  result <- harvest(df, tgt, method="lbfgsb", convergence = list(absolute = 1e-6))
  expect_s3_class(result, "data.frame")
  expect_true("weights" %in% names(result))
  diag <- diagnose_weights(df, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-4))
})

test_that("L-BFGS-B max_weight bound respected", {
  set.seed(1)
  n <- 10000L
  x <- sample(c("a","b"), n, replace=TRUE, prob=c(0.9,0.1))
  df <- data.frame(x=factor(x))
  tgt <- list(x=c(a=0.5, b=0.5))
  result <- harvest(df, tgt, method="lbfgsb", max_weight=5, convergence = list(absolute = 1e-6))
  expect_true(max(result$weights) <= 5.0 + 1e-10)
})

test_that("near-one max_weight rejected for lbfgsb, accepted for ieppa", {
  set.seed(1)
  df <- data.frame(x = factor(sample(c("a","b"), 200, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))

  # lbfgsb: exact and near-1 max_weight → logit singularity error
  expect_error(harvest(df, tgt, method="lbfgsb", max_weight=1.0, convergence = list(absolute = 1e-6)),
               regexp="logit link undefined")
  expect_error(harvest(df, tgt, method="lbfgsb", max_weight=1.0 + 5e-7, convergence = list(absolute = 1e-6)),
               regexp="logit link undefined")
  expect_error(harvest(df, tgt, method="lbfgsb", max_weight=1.0 - 5e-7, convergence = list(absolute = 1e-6)),
               regexp="logit link undefined")

  # ieppa: near-1 max_weight → valid (no logit link); may or may not converge
  expect_no_error(suppressWarnings(
    harvest(df, tgt, method="ieppa", max_weight=1.0 + 5e-7, convergence = list(absolute = 1e-6))
  ))

  # lbfgsb: max_weight=2.0 (well outside eps) → valid
  # suppressWarnings: small n=200 balanced sample may emit convergence warning; not the assertion under test
  expect_no_error(suppressWarnings(harvest(df, tgt, method="lbfgsb", max_weight=2.0, convergence = list(absolute = 1e-6))))
})

test_that("L-BFGS-B converges with tight bounds (max=1.5, min=0.2)", {
  set.seed(11L)
  n   <- 500L
  df  <- data.frame(x = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.3, c = 0.2))
  res <- leafblower::harvest(df, tgt, method = "lbfgsb",
                              max_weight = 1.5, min_weight = 0.2,
                              convergence = list(absolute = 1e-6))
  expect_equal(sum(res$weights), as.double(n), tolerance = 1e-6)
  expect_true(max(res$weights) <= 1.5 + 1e-6)
  expect_true(min(res$weights) >= 0.2 - 1e-6)
})

test_that("L-BFGS-B stable near infeasibility boundary (90/10 split, tight bounds)", {
  set.seed(99L)
  n   <- 300L
  df  <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE,
                                prob = c(0.9, 0.1))))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  res <- leafblower::harvest(df, tgt, method = "lbfgsb",
                              max_weight = 3.0, min_weight = 0.1,
                              convergence = list(absolute = 1e-6))
  expect_equal(sum(res$weights), as.double(n), tolerance = 1e-5)
  expect_true(max(res$weights) <= 3.0 + 1e-5)
  expect_true(min(res$weights) >= 0.1 - 1e-5)
})

test_that("L-BFGS-B sum=n at solver exit, loose bounds", {
  set.seed(12L)
  n   <- 1000L
  df  <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  res <- leafblower::harvest(df, tgt, method = "lbfgsb",
                              max_weight = 100, min_weight = 0,
                              convergence = list(absolute = 1e-6))
  expect_equal(sum(res$weights), as.double(n), tolerance = 1e-6)
})

test_that("B2: lbfgsb emits BUDGET(4) when max_iterations exhausted", {
  set.seed(42)
  n  <- 5000L
  df <- data.frame(x = factor(sample(c("a","b"), n, replace=TRUE, prob=c(0.9,0.1))))
  tgt <- list(x = c(a=0.5, b=0.5))
  result <- suppressWarnings(
    harvest(df, tgt, method="lbfgsb", max_iterations=3L, max_weight=10,
            convergence=list(pct=1e-4))
  )
  s <- attr(result, "result")$status
  expect_equal(s, 4L, info=paste("expected BUDGET=4, got:", s))
})

test_that("B2: lbfgsb never emits NOCONV(1)", {
  set.seed(7)
  n   <- 3000L
  df  <- data.frame(x = factor(sample(c("a","b","c"), n, replace=TRUE)))
  tgt <- list(x = c(a=1/3, b=1/3, c=1/3))
  result <- suppressWarnings(
    harvest(df, tgt, method="lbfgsb", max_iterations=5L)
  )
  s <- attr(result, "result")$status
  expect_false(s == 1L, info=paste("got NOCONV(1), expected BUDGET(4) or OK(0)"))
})

test_that("R5: lbfgsb convergence_used$rule reflects requested rule", {
  set.seed(42)
  n   <- 5000L
  df  <- data.frame(x = factor(sample(c("a","b","c"), n, replace=TRUE)))
  tgt <- list(x = c(a=1/3, b=1/3, c=1/3))
  result <- suppressWarnings(
    harvest(df, tgt, method="lbfgsb", max_weight=5,
            convergence=list(rule="improvement", pct=1e-4))
  )
  cu <- attr(result, "result")$convergence_used
  expect_equal(cu$rule, "improvement",
               info=paste("expected 'improvement', got:", cu$rule))
})
