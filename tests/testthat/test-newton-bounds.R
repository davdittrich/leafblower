# Tests for leafblower-73d7: newton_kl bounds-violation downgrade-first fix.
#
# Fix: any-violation (n_violated > 0) now sets RK_ERR_NOCONV BEFORE writing
# st.weights, so caller never sees violated weights paired with status RK_OK.
#
# DONE_WITH_CONCERNS — fixture limitation:
#   The bounds-violation weight-recovery path (step 6 of newton_calibrate) is
#   only reachable when the Newton inner loop converges (status=RK_OK at end of
#   inner loop) AND the recovered weights w[i] exceed [min_weight, max_weight].
#   This requires max_weight to lie strictly between the last-accepted iterate
#   weight and the final unconstrained weight — a narrow window.
#   Empirical sweep (set.seed(77L), 5-group, max_weight 2.6310..2.6320) shows
#   the inner loop ITSELF diverges (status=1, n_viol=0) for max_weight below
#   the natural solution weight (2.6316), because the dual landscape degenerates
#   when the bound is tighter than the unconstrained solution. No seed/skew
#   combination found that deterministically places 1-5% of obs in the violation
#   window while keeping the inner loop convergent.
#
#   The C fix is correct and verified by code inspection:
#     - src/newton_calib.cpp: `n_violated > 0` replaces `frac_violated > 0.05`
#     - early return skips the st.weights and res.base.best_weights writes
#     - max_error left at default 0.0 (finite) on the violation path
#   The behavioral change is: a 1-5% violation that previously yielded
#   status=RK_OK + violating weights now yields status=RK_ERR_NOCONV + safe
#   pre-violation weights. A follow-up ticket should hand-craft a C-level unit
#   test that injects a controlled violation directly into newton_calibrate().

# T1: Post-fix: the >5% violation path (triggered by inner-loop NOCONV in
# the extreme-skew fixture) still surfaces status=1 and finite max_error.
# This is the existing T4 scenario, repeated here to guard the fix didn't
# regress the NOCONV exit path.
test_that("newton-bounds T1: >5%-violation scenario: status=RK_ERR_NOCONV (=1), finite max_error", {
  set.seed(42L); n <- 2000L
  df  <- data.frame(grp = factor(sample(c("A", "B"), n, TRUE, prob = c(0.5, 0.5))))
  tgt <- list(grp = c(A = 0.95, B = 0.05))
  r <- suppressWarnings(
    harvest(df, tgt, method = "newton_kl", max_weight = 1.3, min_weight = 0,
            max_iterations = 50L, attach_weights = FALSE, verbose = 0L))
  res <- attr(r, "result")
  expect_equal(res$status, 1L,
    label = sprintf("newton-bounds T1: status=%d (expected RK_ERR_NOCONV=1)", res$status))
  expect_true(is.finite(res$max_error),
    label = "newton-bounds T1: max_error must be finite on NOCONV exit")
})

# T2: NOCONV from inner loop divergence: n_bounds_violated == 0 and status != 0.
# The step-6 bounds-violation path is only reached when the inner Newton loop
# converges AND the recovered weights exceed [min_weight, max_weight]. When the
# inner loop diverges first (status=NOCONV before reaching step 6), n_bounds_violated
# remains 0. This distinguishes the two NOCONV exit paths and confirms the fix
# did not accidentally set n_bounds_violated on inner-loop failures.
test_that("newton-bounds T2: inner-loop NOCONV has n_bounds_violated=0 and status!=0", {
  set.seed(42L); n <- 2000L
  df  <- data.frame(grp = factor(sample(c("A", "B"), n, TRUE, prob = c(0.5, 0.5))))
  tgt <- list(grp = c(A = 0.95, B = 0.05))
  r <- suppressWarnings(
    harvest(df, tgt, method = "newton_kl", max_weight = 1.3, min_weight = 0,
            max_iterations = 50L, attach_weights = FALSE, verbose = 0L))
  res <- attr(r, "result")
  # Inner loop diverges before step-6 weight recovery; status != 0.
  expect_true(res$status != 0L,
    label = sprintf("newton-bounds T2: expected non-zero status, got %d", res$status))
  # n_bounds_violated == 0 because weight recovery (step 6) was never reached.
  expect_equal(res$n_bounds_violated, 0L,
    label = sprintf("newton-bounds T2: expected n_bounds_violated=0 (inner-loop exit), got %d",
                    res$n_bounds_violated))
})

# T3: Successful convergence within bounds: status=RK_OK and weights within bounds.
# Verifies the fix didn't break the happy path (n_violated == 0 still writes
# st.weights and res.base.best_weights normally).
test_that("newton-bounds T3: convergent unconstrained run: status=RK_OK and weights within bounds", {
  set.seed(55L); n <- 2000L
  df <- data.frame(
    grp = factor(sample(c("A", "B", "C"), n, TRUE, prob = c(0.60, 0.37, 0.03)))
  )
  tgt <- list(grp = c(A = 0.50, B = 0.47, C = 0.03))
  # max_weight generous enough that Newton converges without violation.
  r <- suppressWarnings(
    harvest(df, tgt, method = "newton_kl", max_weight = 5.0, min_weight = 0,
            max_iterations = 200L, attach_weights = FALSE, verbose = 0L))
  res <- attr(r, "result")
  expect_equal(res$status, 0L,
    label = sprintf("newton-bounds T3: expected RK_OK=0, got %d", res$status))
  expect_true(all(r >= 0 - 1e-12),
    label = "newton-bounds T3: weights must be >= min_weight=0")
  expect_true(all(r <= 5.0 + 1e-12),
    label = sprintf("newton-bounds T3: max weight %.6f > 5.0", max(r)))
  expect_lt(res$max_error, 1e-6,
    label = sprintf("newton-bounds T3: max_error=%.2e must be < 1e-6", res$max_error))
})
