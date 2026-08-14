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

# ---------------------------------------------------------------------------
# Task 2: 50 fixed, stratified, heavy-tailed datasets (KPI-02, D-02, D-03,
# D-04, D-05). Generalizes the tracer above into the full property sweep.
#
# Stratification (D-04) -- 50 datasets split into three groups so a bound
# violation stays attributable to one axis while the interaction case is
# still probed:
#   - "skewed-weights-only" (17, seeds 2001-2017): near-balanced category
#     marginals (no sparse cell), strong lognormal + Cauchy-contaminant
#     design-weight skew drives the clamp.
#   - "sparse-cells-only"   (16, seeds 3001-3016): severely skewed category
#     marginals (a 3% minority cell), mild design-weight skew (small
#     lognormal sigma, light contamination) so sparsity -- not weight
#     skew -- is the axis under stress.
#   - "both"                (17, seeds 4001-4017): severe on both axes; this
#     is the Task 1 tracer's fixture shape, replicated across 17 seeds.
#   17 + 16 + 17 = 50.
#
# Design weights (D-03) are NEVER Gaussian in any stratum: a lognormal bulk
# plus a Cauchy heavy-tailed contaminant on a small row fraction. Gaussian
# design weights never push a calibrated result against [min_weight,
# max_weight], so a bound assertion on Gaussian-drawn weights would be
# vacuous. Parameters were derived to actually reach the clamps, not tuned
# post hoc to pass: measured during planning, ALL 50 datasets engaged the
# clamp on at least one of the 8 swept solvers (50/50), and neutralizing the
# contaminant (a constant x1.0 multiplier in place of the Cauchy draw) drops
# that to 34/50 -- below the 40-dataset floor asserted below, which is what
# proves the witness actually detects a softened generator rather than
# passing regardless.
.bp_gen <- function(n, x_probs, y_probs, sdlog, contam_frac, contam_loc, contam_scale) {
  x <- factor(sample(names(x_probs), n, replace = TRUE, prob = x_probs))
  y <- factor(sample(names(y_probs), n, replace = TRUE, prob = y_probs))
  dw <- rlnorm(n, meanlog = 0, sdlog = sdlog)
  k <- floor(contam_frac * n)
  if (k > 0) dw[sample(n, k)] <- abs(rcauchy(k, location = contam_loc, scale = contam_scale))
  list(df = data.frame(x = x, y = y), dw = dw, target = list(x = x_probs, y = y_probs))
}

# The eight bounds-enforcing solvers (D-01). newton_kl is deliberately NOT
# included -- it has a different shipped reporting contract, adjudicated
# separately in Task 3/Task 4 below.
.bp_solvers <- c("oris", "raking", "sinkhorn", "chebyshev", "greg", "oris_soft",
                  "greenkhorn", "logit")
.bp_min_weight <- 0.2
.bp_max_weight <- 5

.bp_datasets <- list()

# Stratum: skewed-weights-only (17 datasets)
set.seed(2001L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2002L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2003L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2004L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2005L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2006L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2007L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2008L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2009L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2010L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2011L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2012L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2013L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2014L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2015L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2016L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))
set.seed(2017L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "skewed-weights-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.40, B = 0.35, C = 0.25), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.3, contam_frac = 0.06, contam_loc = 20, contam_scale = 15))

# Stratum: sparse-cells-only (16 datasets)
set.seed(3001L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3002L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3003L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3004L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3005L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3006L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3007L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3008L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3009L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3010L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3011L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3012L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3013L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3014L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3015L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))
set.seed(3016L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "sparse-cells-only", data = .bp_gen(n = 1000L, x_probs = c(A = 0.90, B = 0.07, C = 0.03), y_probs = c(P = 0.5, Q = 0.5), sdlog = 0.4, contam_frac = 0.02, contam_loc = 10, contam_scale = 8))

# Stratum: both (17 datasets)
set.seed(4001L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4002L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4003L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4004L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4005L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4006L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4007L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4008L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4009L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4010L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4011L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4012L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4013L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4014L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4015L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4016L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))
set.seed(4017L); .bp_datasets[[length(.bp_datasets) + 1]] <- list(stratum = "both", data = .bp_gen(n = 1200L, x_probs = c(A = 0.80, B = 0.15, C = 0.05), y_probs = c(P = 0.55, Q = 0.45), sdlog = 1.2, contam_frac = 0.05, contam_loc = 20, contam_scale = 15))

stopifnot(length(.bp_datasets) == 50L)

test_that("bound invariant holds across 50 fixed stratified datasets and 8 solvers (KPI-02)", {
  n_engaged <- 0L
  for (ds in .bp_datasets) {
    f <- ds$data
    engaged <- FALSE
    for (m in .bp_solvers) {
      res <- suppressWarnings(
        harvest(f$df, f$target, method = m, design_weights = f$dw,
                min_weight = .bp_min_weight, max_weight = .bp_max_weight,
                bounds_mode = "unit", max_iterations = 1000L, attach_weights = TRUE)
      )
      w <- res$weights
      r <- attr(res, "result")

      # KPI-02, unconditional per D-05: no convergence precheck gates these.
      expect_true(all(w >= .bp_min_weight - 1e-10),
                  info = sprintf("[%s/%s] min(w)=%.15f < %.1f - 1e-10",
                                  ds$stratum, m, min(w), .bp_min_weight))
      expect_true(all(w <= .bp_max_weight + 1e-10),
                  info = sprintf("[%s/%s] max(w)=%.15f > %.1f + 1e-10",
                                  ds$stratum, m, max(w), .bp_max_weight))

      if (!is.null(r$n_bounds_clamped) && r$n_bounds_clamped > 0L) engaged <- TRUE
    }
    if (engaged) n_engaged <- n_engaged + 1L
  }

  # Non-vacuity witness: a dataset counts as "engaged" if the clamp fired on
  # at least one of the 8 swept solvers. Measured during planning at 50/50;
  # a floor of 40 has headroom while still failing loudly if a future edit
  # softens the generator (verified: neutralizing the contaminant drops the
  # count to 34/50, below this floor).
  expect_gte(n_engaged, 40L)
})
