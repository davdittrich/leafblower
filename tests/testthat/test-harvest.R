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

test_that("convergence$absolute is forwarded to solver", {
  set.seed(42)
  n <- 10000L
  df <- data.frame(
    age = factor(sample(c("Y","M","O"), n, replace=TRUE,
                        prob=c(0.60, 0.30, 0.10))),
    sex = factor(sample(c("M","F"),     n, replace=TRUE,
                        prob=c(0.70, 0.30)))
  )
  tgt <- list(
    age = c(Y=0.33, M=0.40, O=0.27),
    sex = c(M=0.49, F=0.51)
  )
  # 2 iterations with default tol (1e-6): competing margins cannot converge -> warning
  expect_warning(
    harvest(df, tgt, method="ieppa", max_iterations=2),
    regexp="did not converge"
  )
  # 2 iterations with loose tol (0.3): error after 2 iters < 0.3 -> no warning
  # Before fix: tol_abs ignored, tol=1e-6 used, warning fires -> test fails
  # After fix:  tol_abs=0.3 forwarded, error < 0.3 accepted -> no warning
  expect_no_warning(
    harvest(df, tgt, method="ieppa", max_iterations=2,
            convergence=list(absolute=0.3))
  )
})
