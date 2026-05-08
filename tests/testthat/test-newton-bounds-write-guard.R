# test-newton-bounds-write-guard.R
#
# Verifies the write-guard introduced in leafblower-73d7:
# when frac_violated > 0.05 (T4 contract), newton_kl must NOT write the
# violating weight vector to res.base.best_weights.  The C++ guard leaves
# best_weights empty; r_bridge fills it with sentinel zeros (length n).
# The caller therefore observes status = RK_ERR_NOCONV AND a weight vector
# that is entirely zero — not the out-of-bounds iterate.
#
# Fixture: same extreme-skew setup as T4 (50/50 sample, 95/5 target,
# max_weight = 1.3).  Unconstrained Newton weights for "A" ~ 1.9, so ~50 %
# exceed the bound — well above the 5 % gate.

test_that("write-guard: status=RK_ERR_NOCONV and best_weights is sentinel zeros on >5% violation", {
  set.seed(42L); n <- 2000L
  df  <- data.frame(grp = factor(sample(c("A","B"), n, TRUE, prob=c(0.5, 0.5))))
  tgt <- list(grp = c(A = 0.95, B = 0.05))

  r   <- suppressWarnings(
    harvest(df, tgt, method = "newton_kl",
            max_weight = 1.3, min_weight = 0,
            max_iterations = 50L, attach_weights = FALSE, verbose = 0L)
  )
  res <- attr(r, "result")

  # T4 contract: status == RK_ERR_NOCONV
  expect_equal(res$status, 1L,
    label = sprintf("write-guard: status=%d (expected RK_ERR_NOCONV=1)", res$status))

  # n_bounds_violated is now surfaced on the direct newton_kl path
  expect_true(res$n_bounds_violated > 0L,
    label = "write-guard: n_bounds_violated must be > 0 when threshold crossed")

  # Write-guard: best_weights must be the sentinel zero vector (all zeros),
  # NOT the violating iterate.  The C++ guard leaves best_weights empty;
  # r_bridge fills with zeros (length n) as the safe sentinel.
  bw <- res$best_weights
  expect_equal(length(bw), n,
    label = sprintf("write-guard: best_weights length should be n=%d", n))
  expect_true(all(bw == 0),
    label = "write-guard: best_weights must be all-zero sentinel when violation guard triggers")
})
