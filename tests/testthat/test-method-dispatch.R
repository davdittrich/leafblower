library(leafblower)

methods_to_test <- c("oris", "oris_soft", "raking", "greg",
                     "chebyshev", "sinkhorn", "auto", "greenkhorn", "logit")

for (m in methods_to_test) {
  local({
    method <- m
    test_that(paste("method dispatch works for", method), {
      n   <- 20L
      df  <- data.frame(x = factor(rep(c("a", "b"), n / 2L)))
      tgt <- list(x = c(a = 0.5, b = 0.5))
      expect_no_error(
        harvest(df, tgt, method = method, max_iterations = 10L,
                convergence = list(absolute = 0.1))
      )
    })
  })
}

test_that("lbfgsb method is rejected with match.arg error", {
  df  <- data.frame(a = factor(c("x","y","x"), levels = c("x","y")))
  tgt <- list(a = c(x = 0.6, y = 0.4))
  expect_error(harvest(df, tgt, method = "lbfgsb"))
})

test_that("AUTO path completes quickly on hard input (no lbfgsb hang)", {
  # Severe skew: data is 95% x but target is 90% y — hard for primary solver
  df  <- data.frame(
    a = factor(rep(c("x","y"), c(950, 50)), levels = c("x","y"))
  )
  tgt <- list(a = c(x = 0.1, y = 0.9))
  t0      <- proc.time()["elapsed"]
  r       <- suppressWarnings(
    harvest(df, tgt, method = "auto", max_iterations = 10L, verbose = 0L)
  )
  elapsed <- proc.time()["elapsed"] - t0
  expect_lt(elapsed, 15.0)
  alg     <- attr(r, "result")$algorithm
  expect_true(alg %in% c("oris", "newton_kl", "raking"))
})

test_that("unknown method rejected by R layer (C++ kAlgMap fallback not reachable via harvest)", {
  # harvest.R validates method via match.arg before reaching C++.
  # This test confirms the R-level guard works (catching unknown methods with
  # "should be one of" error). The C++ RK_ALG_RAKING fallback is a safety net
  # for direct C API callers only and is intentionally not tested from R.
  n   <- 20L
  df  <- data.frame(x = factor(rep(c("a", "b"), n / 2L)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  expect_error(
    harvest(df, tgt, method = "nonexistent_method", max_iterations = 10L,
            convergence = list(absolute = 0.1)),
    regexp = "should be one of"
  )
})
