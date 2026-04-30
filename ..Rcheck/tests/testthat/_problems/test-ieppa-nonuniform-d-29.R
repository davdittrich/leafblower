# Extracted from test-ieppa-nonuniform-d.R:29

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "leafblower", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
context("ieppa faithful — non-uniform design weights within cell")

# test -------------------------------------------------------------------------
set.seed(42)
n <- 1000
df <- data.frame(
    a = factor(rep(0:1, each = n/2)),
    b = factor(rep(0L, n))
  )
d <- rep(c(1, 10), length.out = n)
tgt <- list(
    a = c(`0` = 0.3, `1` = 0.7),
    b = c(`0` = 1.0)
  )
res <- harvest(df, tgt, method = "ieppa",
                 start_weights = d, max_weight = 2.0, min_weight = 0,
                 convergence = list(absolute = 1e-6))
w <- res$weights
m_a0 <- sum(w[df$a == 0]) / sum(w)
m_a1 <- sum(w[df$a == 1]) / sum(w)
expect_lt(abs(m_a0 - 0.3), 1e-4)
expect_lt(abs(m_a1 - 0.7), 1e-4)
