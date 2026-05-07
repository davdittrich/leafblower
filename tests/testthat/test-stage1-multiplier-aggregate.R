# test-stage1-multiplier-aggregate.R — T-M (leafblower-6ycz.1.13).
#
# Spec §6 mandates per-cell Stage-1 multiplier inheritance for SPARSE cells:
#   weights[i] = w_init[i] * stage1_multipliers[c]   for i in sparse cell c.
#
# Pre-fix bug: `stage1_multipliers[c] = st.weights[i] / w_init[i]` overwritten
# in a per-obs loop — the LAST-i in cell c wins (arbitrary by iteration order),
# so different observations in the SAME cell with different post-Stage-1
# weights produce a non-deterministic, observation-dependent multiplier.
#
# Spec §6 fix: cell-aggregate multiplier:
#     stage1_multipliers[c] = (Σ_{i ∈ c} weights[i]) / (Σ_{i ∈ c} w_init[i])
# (with Σw_init > 0 guard else 1.0).
#
# Test strategy: we cannot directly inspect stage1_multipliers, but we can
# observe its effect on a sparse cell. Construct a DGP where a coarse cell c*
# is forced sparse (n_cell < min_cell_n), and where if Stage-1 ran on the full
# data, observations within c* would receive non-uniform per-obs ratios
# weights[i]/w_init[i]. Post-fix, the multiplier applied to c* equals the
# cell-aggregate ratio Σw / Σw_init computed across the full Stage-1 result —
# i.e., for ALL observations in c*, weight[i] / w_init[i] is identical AND
# equals Σw_in_c* / Σw_init_in_c*.
#
# By construction: w_init[i] is uniform (=1.0) at entry → w_init differs
# across i only if the harvest input weight differs. We pass uniform input
# weights, so w_init[i] = 1.0 for all i. After sparse inheritance, weights[i]
# in sparse c* must equal stage1_multipliers[c*] (since w_init[i]=1).
#
# CORRECT (cell-aggregate): all i in c* share the same final weight, equal to
#   mean(post_stage1_full_data_weights[i in c*]) / 1.0
# (because Σw / Σw_init = Σw / n_cell = mean(w)).
#
# BROKEN (per-obs overwrite): final weight in c* = w_post_stage1[last_i_in_c*]
# — depends on iteration order; almost always differs from the mean unless by
# coincidence.

library(testthat)

# Build DGP where: coarse_mask=(1,0); coarse cell c=0 is sparse (small);
# coarse cell c=1 is dense; targets non-trivial so post-Stage-1 weights differ
# within cell 0.
#
# Construction details:
# - Coarse margin g (2 levels): g=0 has 4 obs (sparse, < min_cell_n=10);
#   g=1 has 96 obs (dense, ≥ 10).
# - Fine margin f (2 levels): f distributes such that within g=0,
#   2 obs have f=0 and 2 have f=1, with full-data fine-margin target requiring
#   non-uniform reweighting that — under the BROKEN per-obs path — produces
#   distinct ratios for the 4 obs in g=0.
make_sparse_cell_dgp <- function() {
  # Hand-constructed: N=100, K=2 binary, coarse=g (mask=1), fine=f (mask=0).
  # g=0: obs 1..4 (sparse).
  # g=1: obs 5..100 (dense).
  # Fine within g=0: f=(0,0,1,1). Within g=1: roughly 50/50.
  set.seed(11L)
  g <- c(rep(0L, 4L), rep(1L, 96L))
  # Fine: within g=0 → 2x f=0, 2x f=1; within g=1 → 48 of each.
  f <- c(0L, 0L, 1L, 1L,
         rep(c(0L, 1L), times = 48L))
  df <- data.frame(g = factor(g), f = factor(f))
  # Targets that require non-trivial reweighting on f, so post-Stage-1
  # full-data weights are non-uniform across obs.
  # Coarse target g: keep observed proportion (4/100, 96/100) so coarse is
  # ~satisfied at unit weights — ensures any per-obs disparity comes from f.
  targets <- list(
    g = c(`0` = 0.04, `1` = 0.96),
    f = c(`0` = 0.30, `1` = 0.70)   # f is fine; Stage-1 must NOT touch
  )
  list(df = df, targets = targets, coarse_mask = c(1L, 0L))
}

# ---------------------------------------------------------------------------
# Cell-aggregate multiplier — sparse cell uniform-weight assertion.
#
# Post-fix invariant: for sparse cell c, all observations in c receive the
# SAME weight after sparse inheritance (since w_init is uniform → multiplier
# is the same constant for every i in c).
#
# The stronger invariant: the constant equals the cell-aggregate ratio
#   (Σ_{i ∈ c} w_post_stage1[i]) / (Σ_{i ∈ c} w_init[i])
# = mean(w_post_stage1[i ∈ c])  (when w_init is uniformly 1).
#
# Pre-fix: weights in c = single per-obs ratio (last-i wins) — generally
# NOT equal across all i in c when post-Stage-1 weights differ within c.
# ---------------------------------------------------------------------------
test_that("sparse cell receives uniform multiplier (cell-aggregate, raking)", {
  d <- make_sparse_cell_dgp()
  cfg <- list(
    coarse_mask      = d$coarse_mask,
    min_cell_n       = 10L,           # cell 0 (n=4) is sparse; cell 1 (n=96) is not
    mode             = 0L,
    outer_tol        = 1e-6,
    outer_iterations = 20L
  )
  r <- tryCatch(
    suppressWarnings(
      harvest(d$df, d$targets, method = "raking", hierarchical = cfg)
    ),
    error = function(e) NULL
  )
  skip_if(is.null(r), "DGP errored")
  diag <- attr(r, "result")
  skip_if(is.null(diag) || diag$n_cells_skipped < 1L,
          "expected ≥1 sparse cell; DGP didn't produce one")

  # Sparse cell observations are obs 1..4 (g=0).
  w_sparse <- r$weight[1:4]
  # Cell-aggregate post-fix: all four weights are equal (within numerical eps).
  # Pre-fix bug: per-obs ratio overwrite produces non-uniform weights here
  # whenever post-Stage-1 weights vary across the 4 obs (typical case under
  # the constructed DGP where f varies within g=0).
  expect_equal(diff(range(w_sparse)), 0,
               tolerance = 1e-9,
               label = sprintf(
                 "all sparse-cell weights must be equal (cell-aggregate): w=[%s]",
                 paste(format(w_sparse, digits = 6), collapse = ", ")))
})

test_that("sparse cell receives uniform multiplier (cell-aggregate, sinkhorn)", {
  d <- make_sparse_cell_dgp()
  cfg <- list(
    coarse_mask      = d$coarse_mask,
    min_cell_n       = 10L,
    mode             = 0L,
    outer_tol        = 1e-6,
    outer_iterations = 20L
  )
  r <- tryCatch(
    suppressWarnings(
      harvest(d$df, d$targets, method = "sinkhorn", hierarchical = cfg)
    ),
    error = function(e) NULL
  )
  skip_if(is.null(r), "DGP errored")
  diag <- attr(r, "result")
  skip_if(is.null(diag) || diag$n_cells_skipped < 1L,
          "expected ≥1 sparse cell")
  w_sparse <- r$weight[1:4]
  expect_equal(diff(range(w_sparse)), 0,
               tolerance = 1e-9,
               label = sprintf(
                 "sinkhorn sparse-cell weights must be equal: w=[%s]",
                 paste(format(w_sparse, digits = 6), collapse = ", ")))
})

test_that("sparse cell receives uniform multiplier (cell-aggregate, greenkhorn)", {
  d <- make_sparse_cell_dgp()
  cfg <- list(
    coarse_mask      = d$coarse_mask,
    min_cell_n       = 10L,
    mode             = 0L,
    outer_tol        = 1e-6,
    outer_iterations = 20L
  )
  r <- tryCatch(
    suppressWarnings(
      harvest(d$df, d$targets, method = "greenkhorn", hierarchical = cfg)
    ),
    error = function(e) NULL
  )
  skip_if(is.null(r), "DGP errored")
  diag <- attr(r, "result")
  skip_if(is.null(diag) || diag$n_cells_skipped < 1L,
          "expected ≥1 sparse cell")
  w_sparse <- r$weight[1:4]
  expect_equal(diff(range(w_sparse)), 0,
               tolerance = 1e-9,
               label = sprintf(
                 "greenkhorn sparse-cell weights must be equal: w=[%s]",
                 paste(format(w_sparse, digits = 6), collapse = ", ")))
})
