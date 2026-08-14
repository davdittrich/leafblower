library(leafblower)

# KPI-02 / ROADMAP SC4: `min(w) >= min_weight` and `max(w) <= max_weight`
# within 1e-10, asserted unconditionally on the returned weight vector
# regardless of solver status (D-05). Bounds are enforced in the shared
# C++ core (`lbw::finalize_weights_buf`), so R testthat alone exercises
# the same code path both bindings share (D-01).
#
# `bounds_mode = "unit"` is mandatory, not a style choice: the default
# cell-aggregate mode only COUNTS bound violations, it does not enforce
# them (R/harvest.R:833-835 warns exactly this). A property test against
# the default mode would be correctly sometimes-failing, not a defect
# signal (Pattern 2, 01-RESEARCH.md:186-206).
#
# This task proves the call shape, result-field access, and assertion on
# ONE dataset / ONE solver end to end. Task 2 generalizes this into a
# fixed 50-dataset stratified sweep across the eight bounds-enforcing
# solvers.

# A skewed-category / heavy-tailed-design-weight fixture measured during
# planning: n=1200, a three-level factor at 0.80/0.15/0.05, a two-level
# factor at 0.55/0.45, lognormal design weights at sdlog=1.2, a Cauchy
# contaminant on ~5% of rows, bounds [0.2, 5]. Under raking that fixture
# drives BOTH clamps to bind exactly (min(w) == 0.2, max(w) == 5).
.bound_fixture_1 <- function() {
  set.seed(7L)
  n <- 1200L
  x <- factor(sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.80, 0.15, 0.05)))
  y <- factor(sample(c("P", "Q"), n, replace = TRUE, prob = c(0.55, 0.45)))
  # Lognormal bulk (D-03): mostly well-behaved design weights.
  dw <- rlnorm(n, meanlog = 0, sdlog = 1.2)
  # Heavy-tailed contaminant fraction (D-03): a Cauchy overlay on ~5% of
  # rows, never a Gaussian draw. Gaussian design weights never push the
  # calibrated result against the clamps, so a bound assertion on such
  # data would be vacuous (D-03).
  k <- floor(0.05 * n)
  dw[sample(n, k)] <- abs(rcauchy(k, location = 20, scale = 15))
  list(
    df = data.frame(x = x, y = y),
    dw = dw,
    target = list(x = c(A = 0.80, B = 0.15, C = 0.05), y = c(P = 0.55, Q = 0.45))
  )
}

test_that("bound invariant holds end to end: one dataset, one solver (KPI-02 tracer)", {
  f <- .bound_fixture_1()
  res <- suppressWarnings(
    harvest(f$df, f$target, method = "raking", design_weights = f$dw,
            min_weight = 0.2, max_weight = 5, bounds_mode = "unit",
            max_iterations = 1000L, attach_weights = TRUE)
  )
  w <- res$weights
  r <- attr(res, "result")

  # KPI-02, literal wording, asserted unconditionally regardless of status
  # (D-05): even a non-converged iterate must respect its bounds.
  expect_true(all(w >= 0.2 - 1e-10),
              info = sprintf("min(w)=%.15f < 0.2 - 1e-10", min(w)))
  expect_true(all(w <= 5 + 1e-10),
              info = sprintf("max(w)=%.15f > 5 + 1e-10", max(w)))

  # Non-vacuity: the clamp counter on the result attribute must be > 0, so
  # the two assertions above are proven to actually exercise the water-fill
  # rather than passing on already-feasible weights.
  expect_gt(r$n_bounds_clamped, 0L)
})
