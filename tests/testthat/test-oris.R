context("oris (faithful algBCD)")

test_that("ORIS converges: 1 margin, 2 cats, no bounds", {
  set.seed(42)
  n   <- 100L
  df  <- data.frame(x = factor(sample(c("a","b"), n, replace=TRUE, prob=c(0.7,0.3))))
  tgt <- list(x = c(a=0.5, b=0.5))
  result <- harvest(df, tgt, method="oris", convergence = list(absolute = 1e-6))
  expect_true(attr(result, "algorithm") == "oris")
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("ORIS respects max_weight=2 on tight bounds", {
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
  result <- harvest(df, tgt, method="oris", max_weight=2, convergence = list(absolute = 1e-6))
  expect_true(max(result$weights) <= 2.0 + 1e-8)
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("ORIS respects min_weight=0.5", {
  set.seed(3)
  n   <- 10000L
  df  <- data.frame(
    x = factor(sample(c("a","b","c","d","e"), n, replace=TRUE))
  )
  tgt <- list(x = c(a=0.2, b=0.2, c=0.2, d=0.2, e=0.2))
  result <- harvest(df, tgt, method="oris", min_weight=0.5, max_weight=5, convergence = list(absolute = 1e-6))
  expect_true(min(result$weights) >= 0.5 - 1e-8)
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("ORIS output weights have mean=1 and respect bounds", {
  set.seed(5L)
  n   <- 1000L
  df  <- data.frame(x = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.3, c = 0.2))
  res <- leafblower::harvest(df, tgt, method = "oris",
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
    method = "oris",
    max_iterations = 5L,
    convergence = list(rule="improvement", pct=0.5)
  )
  # With correct X_prev initialization, already-calibrated input converges at iter=1.
  expect_lte(attr(result,"result")$iterations, 3L)
})

test_that("B12: oris greedy scheduler produces finite errRp (not 0 sentinel) on non-trivial input", {
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
    method = "oris",
    scheduler = "greedy",
    max_iterations = 50L,
    convergence = list(absolute = 1e-4)
  )
  # Greedy scheduler must converge to correct solution, not silently return 0.
  expect_lt(attr(result,"result")$max_error, 1e-4)
})

# ──────────────────────────────────────────────────────────────────────────────
# CR-A6 (mxcl.6): the SRAA (accelerate=TRUE) path must report max_error consistent
# with the RETURNED, capacity-clamped weights. The SRAA sweep is unconstrained
# (capacity enforced only by the level-exit mass-preserving clamp), so
# res.base.max_error = r.err_rp (oris.cpp:877) is the UNCONSTRAINED pre-clamp
# errRp. On a tight/infeasible problem it reported ~2e-16 ("converged") while the
# returned weights miss margins by ~0.26 — a silent false-convergence. The
# non-SRAA (accelerate=FALSE) path already reports the honest constrained residual.
# ──────────────────────────────────────────────────────────────────────────────

test_that("oris SRAA reports max_error consistent with returned weights (CR-A6)", {
  set.seed(99)
  n <- 1000L
  # n/M_cell ~ 83 -> log SRAA path under natural dispatch; tight max_weight makes
  # the unconstrained sweep and the bounds-clamped output diverge.
  x <- factor(sample(c("H", "L"), n, replace = TRUE, prob = c(0.7, 0.3)))
  y <- factor(sample(c("P", "Q"), n, replace = TRUE, prob = c(0.3, 0.7)))
  z <- factor(sample(c("M", "N", "O"), n, replace = TRUE, prob = c(0.5, 0.3, 0.2)))
  data <- data.frame(x = x, y = y, z = z)
  target <- list(x = c(H = 0.4, L = 0.6),
                 y = c(P = 0.7, Q = 0.3),
                 z = c(M = 0.4, N = 0.35, O = 0.25))

  w <- suppressWarnings(harvest(data, target, method = "oris", accelerate = TRUE,
                                max_iterations = 500L, min_weight = 0.5, max_weight = 1.5,
                                attach_weights = FALSE, verbose = 0))
  r <- attr(w, "result")

  # Achieved max margin error recomputed from the RETURNED weights.
  Wt <- sum(w)
  achieved <- max(vapply(names(target), function(v) {
    max(abs(tapply(w, data[[v]], sum)[names(target[[v]])] / Wt - target[[v]]))
  }, numeric(1)))

  # Reported max_error must reflect the returned weights, not the unconstrained
  # pre-clamp iterate (bug reported ~2e-16 while achieved ~0.26).
  expect_lt(abs(r$max_error - achieved), 1e-6)
})
