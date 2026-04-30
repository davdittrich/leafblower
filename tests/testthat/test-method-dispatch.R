library(leafblower)

methods_to_test <- c("ieppa", "ieppa_soft", "lbfgsb", "raking", "greg",
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

test_that("unknown method rejected by R layer (C++ kAlgMap fallback not reachable via harvest)", {
  # harvest.R validates method via match.arg before reaching C++.
  # This confirms the R-level guard works; the C++ fallback (RK_ALG_IEPPA)
  # is a safety net for direct C_rk_calibrate callers only.
  n   <- 20L
  df  <- data.frame(x = factor(rep(c("a", "b"), n / 2L)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  expect_error(
    harvest(df, tgt, method = "nonexistent_method", max_iterations = 10L,
            convergence = list(absolute = 0.1)),
    regexp = "should be one of"
  )
})
