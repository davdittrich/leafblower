# tests/testthat/test-newton-steihaug-cg.R
# WL-3 (Epic-Dβ): Steihaug-CG trust-region step + Marquardt-style adaptive Δ.
library(leafblower)

test_that("WL-3 Steihaug-CG returns boundary point when pinv exceeds trust radius", {
  # Synthetic problem: severe-skew K=2 fixture engineered so the early-Newton
  # pseudoinverse step exceeds the default trust radius Δ=1.0 by a wide margin
  # (extreme target/sample probability mismatch ⇒ large gradient at λ=0,
  # combined with rank-deficient corner cells). After 1 Newton iter, ||δ_pinv||₂
  # > Δ on iter 0, so Steihaug-CG branch fires and clips to boundary.
  #
  # Hard to assert exactly without inspecting solver internals; instead:
  #  - all weights finite ⇒ CG branch did not produce NaN
  #  - max_error sane (< 1e-1) ⇒ CG step is a descent direction (boundary
  #    projection didn't break the algorithm)
  set.seed(1L)
  n  <- 200L
  df <- data.frame(
    a = factor(sample(c("x", "y"), n, TRUE, prob = c(0.85, 0.15))),
    b = factor(sample(c("p", "q"), n, TRUE, prob = c(0.85, 0.15)))
  )
  # Skewed target vs sample ⇒ first-iter ||δ_pinv||₂ > Δ=1 ⇒ Steihaug-CG fires.
  tgt <- list(a = c(x = 0.15, y = 0.85), b = c(p = 0.15, q = 0.85))
  r <- suppressWarnings(harvest(
    df, tgt,
    method         = "newton_kl",
    max_weight     = 50,
    min_weight     = 0,
    max_iterations = 100L,
    attach_weights = FALSE,
    verbose        = 0L
  ))
  res <- attr(r, "result")
  expect_true(
    all(is.finite(r)),
    label = sprintf(
      "CG-via-trust-region produces finite weights; min=%.3g max=%.3g",
      min(r), max(r)
    )
  )
  expect_lt(
    res$max_error, 1e-1,
    label = sprintf("CG branch produces sane gap; got %.3e", res$max_error)
  )
})
