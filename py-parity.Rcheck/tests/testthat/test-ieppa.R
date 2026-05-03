context("ieppa (faithful algBCD)")

test_that("iEPPA converges: 1 margin, 2 cats, no bounds", {
  set.seed(42)
  n   <- 100L
  df  <- data.frame(x = factor(sample(c("a","b"), n, replace=TRUE, prob=c(0.7,0.3))))
  tgt <- list(x = c(a=0.5, b=0.5))
  result <- harvest(df, tgt, method="ieppa", convergence = list(absolute = 1e-6))
  expect_true(attr(result, "algorithm") == "ieppa")
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("iEPPA respects max_weight=2 on tight bounds", {
  set.seed(7)
  n   <- 10000L
  df  <- data.frame(
    age = factor(sample(c("Y","M","O"), n, replace=TRUE, prob=c(0.40,0.33,0.27))),
    sex = factor(sample(c("M","F"), n, replace=TRUE, prob=c(0.55,0.45))),
    edu = factor(sample(c("HS","Col","Grad"), n, replace=TRUE, prob=c(0.4,0.4,0.2)))
  )
  tgt <- list(
    age = c(Y=0.33, M=0.34, O=0.33),
    sex = c(M=0.50, F=0.50),
    edu = c(HS=0.35, Col=0.45, Grad=0.20)
  )
  result <- harvest(df, tgt, method="ieppa", max_weight=2, convergence = list(absolute = 1e-6))
  expect_true(max(result$weights) <= 2.0 + 1e-8)
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("iEPPA respects min_weight=0.5", {
  set.seed(3)
  n   <- 10000L
  df  <- data.frame(
    x = factor(sample(c("a","b","c","d","e"), n, replace=TRUE))
  )
  tgt <- list(x = c(a=0.2, b=0.2, c=0.2, d=0.2, e=0.2))
  result <- harvest(df, tgt, method="ieppa", min_weight=0.5, max_weight=5, convergence = list(absolute = 1e-6))
  expect_true(min(result$weights) >= 0.5 - 1e-8)
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("iEPPA output weights have mean=1 and respect bounds", {
  set.seed(5L)
  n   <- 1000L
  df  <- data.frame(x = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.3, c = 0.2))
  res <- leafblower::harvest(df, tgt, method = "ieppa",
                              max_weight = 2.0, min_weight = 0.2,
                              convergence = list(absolute = 1e-6),
                              attach_weights = FALSE)
  # mean=1 is guaranteed by both the old fixup loop and the new Dykstra projection;
  # this test guards against regressions in the P2 refactor.
  expect_equal(mean(res), 1.0, tolerance = 1e-10)
  expect_true(max(res) <= 2.0 + 1e-8)
  expect_true(min(res) >= 0.2 - 1e-8)
})

test_that("B11: X_prev initialized from X_init not zeros on first homotopy level", {
  # A nearly-converged problem: targets match empirical proportions.
  # With bug: X_prev=zeros → pct_change huge on iter=1 → convergence check corrupted.
  # With fix: X_prev=X_init → pct_change reflects real change from X_init.
  result <- harvest(
    data.frame(x = factor(c("A","A","B","B","B")), w=rep(1,5)),
    target = list(x = c(A=0.4, B=0.6)),
    method = "ieppa",
    max_iterations = 5L,
    convergence = list(rule="improvement", pct=0.5)
  )
  # With correct X_prev initialization, already-calibrated input converges at iter=1.
  expect_lte(attr(result,"result")$iterations, 3L)
})

test_that("B12: ieppa greedy scheduler produces finite errRp (not 0 sentinel) on non-trivial input", {
  # With bug: compute_margin_errRp_linear/log returned 0.0 when W_total<=0,
  # signalling false perfect convergence to greedy scheduler.
  # With fix: returns Inf, so greedy correctly selects max-error margin.
  # Use n=100 with skewed marginals so the problem is non-trivial.
  set.seed(42)
  n <- 100L
  df <- data.frame(
    x = factor(sample(c("A","B"), n, replace=TRUE, prob=c(0.7,0.3))),
    y = factor(sample(c("P","Q"), n, replace=TRUE, prob=c(0.6,0.4)))
  )
  result <- harvest(
    df,
    target = list(x=c(A=0.5,B=0.5), y=c(P=0.5,Q=0.5)),
    method = "ieppa",
    scheduler = "greedy",
    max_iterations = 50L,
    convergence = list(absolute = 1e-4)
  )
  # Greedy scheduler must converge to correct solution, not silently return 0.
  expect_lt(attr(result,"result")$max_error, 1e-4)
})
