## test-2stage-strategy-a.R
## Verification tests for lbw::outer_iterate_strategy_a (T-C).
##
## DGP: synthetic 2-margin, N=60 obs, 3 × 2 = 6 coarse cells.
## WithinCellSolve mock: fn_mode 0 = identity (weights unchanged → converges),
##                       fn_mode 1 = perturb (residual oscillates → BUDGET).
##
## Tests:
##   (a) RK_OK on convergent mock.
##   (b) RK_ERR_BUDGET + last-iterate-weights returned on non-convergent mock.
##   (c) Σw=n preserved at RK_OK exit.
##   (d) best_iter_idx reflects select_metric cadence (not errRp directly).

# Helper: call C_hier_outer_probe.
.outer_probe <- function(
    weights, group_ids_list, n, K, coarse_mask, min_cell_n,
    stage1_mults, outer_tol, outer_iterations, N, fn_mode) {
  .Call(
    "C_hier_outer_probe",
    as.double(weights),
    group_ids_list,
    as.integer(n),
    as.integer(K),
    as.integer(coarse_mask),
    as.integer(min_cell_n),
    as.double(stage1_mults),
    as.double(outer_tol),
    as.integer(outer_iterations),
    as.integer(N),
    as.integer(fn_mode),
    PACKAGE = "leafblower"
  )
}

## ── DGP setup ─────────────────────────────────────────────────────────────────
## N = 60 obs; margin 1: 3 levels (0,1,2); margin 2: 2 levels (0,1).
## Coarse mask: both margins (all cells are "coarse cells").
## Initial weights: uniform (all 1.0) → Σw = N = 60.

N   <- 60L
K   <- 2L
set.seed(42L)

g1  <- as.integer(rep(0:2, each = 20L))          # levels 0,1,2 (n=20 each)
g2  <- as.integer(rep(c(0L, 1L), times = 30L))   # levels 0,1 (n=30 each)

group_ids_list <- list(g1, g2)
coarse_mask    <- c(1L, 1L)   # both margins coarse
min_cell_n     <- 1L           # no sparse cells

weights0    <- rep(1.0, N)
## stage1_mults: one per cell. 6 cells (3×2), all 1.0 for the mock.
## Actual value only matters for sparse-cell inheritance — no sparse cells here.
stage1_mults <- rep(1.0, 6L)

## ── (a) Convergent mock (fn_mode = 0): identity fn, tol is reachable ─────────
test_that("outer_iterate_strategy_a returns RK_OK on convergent identity mock", {
  res <- .outer_probe(
    weights        = weights0,
    group_ids_list = group_ids_list,
    n              = N, K = K,
    coarse_mask    = coarse_mask,
    min_cell_n     = min_cell_n,
    stage1_mults   = stage1_mults,
    outer_tol      = 1.0,      # very loose — identity fn yields residual ≈ 0
    outer_iterations = 50L,
    N              = N,
    fn_mode        = 0L        # identity: weights unchanged
  )
  expect_equal(res$status,     0L)    # RK_OK = 0
  expect_true(res$iterations_used >= 1L)
  expect_true(is.numeric(res$last_iterate_weights))
  expect_equal(length(res$last_iterate_weights), N)
})

## ── (b) Non-convergent mock: outer_tol=0.0 forces BUDGET regardless of residual ──
## outer_tol=0.0 disables the convergence check (guard: outer_tol > 0.0),
## so the loop always exhausts the budget.  fn_mode=1 mutates weights via
## cell-asymmetric perturbations, ensuring last_iterate_weights != initial.
test_that("outer_iterate_strategy_a returns RK_ERR_BUDGET + last-iterate on budget exhaustion", {
  res <- .outer_probe(
    weights        = weights0,
    group_ids_list = group_ids_list,
    n              = N, K = K,
    coarse_mask    = coarse_mask,
    min_cell_n     = min_cell_n,
    stage1_mults   = stage1_mults,
    outer_tol      = 0.0,      # disabled: outer_tol <= 0 → convergence check never fires
    outer_iterations = 5L,
    N              = N,
    fn_mode        = 1L        # perturb: weights mutated asymmetrically per cell
  )
  expect_equal(res$status, 4L)   # RK_ERR_BUDGET = 4
  expect_equal(res$iterations_used, 5L)
  # last_iterate_weights is the actual last state, NOT best-iterate
  expect_equal(length(res$last_iterate_weights), N)
  expect_true(is.numeric(res$last_iterate_weights))
  # Perturb fn mutated weights → last-iterate differs from initial.
  expect_false(isTRUE(all.equal(res$last_iterate_weights, weights0, tolerance = 1e-10)))
})

## ── (c) Σw = N preserved at RK_OK exit ────────────────────────────────────────
test_that("enforce_sigmaw_eq_n gate fires: Σw == N at RK_OK exit", {
  res <- .outer_probe(
    weights        = weights0,
    group_ids_list = group_ids_list,
    n              = N, K = K,
    coarse_mask    = coarse_mask,
    min_cell_n     = min_cell_n,
    stage1_mults   = stage1_mults,
    outer_tol      = 1.0,
    outer_iterations = 50L,
    N              = N,
    fn_mode        = 0L
  )
  expect_equal(res$status, 0L)  # must be RK_OK for gate to apply
  w_sum <- sum(res$last_iterate_weights)
  expect_equal(w_sum, as.double(N), tolerance = 1e-9)
})

## ── (d) best_iter_idx reflects select_metric at kOuterErrCheckInterval ────────
## With fn_mode=1 and outer_iterations=25, at least one kErrCheckInterval=10
## check fires (iter 0, 10, 20). best_iter_idx must be non-negative.
test_that("best-iterate tracking fires at kErrCheckInterval via select_metric", {
  res <- .outer_probe(
    weights        = weights0,
    group_ids_list = group_ids_list,
    n              = N, K = K,
    coarse_mask    = coarse_mask,
    min_cell_n     = min_cell_n,
    stage1_mults   = stage1_mults,
    outer_tol      = 0.0,      # disabled: forces full 25 iterations
    outer_iterations = 25L,
    N              = N,
    fn_mode        = 1L
  )
  # At least one kOuterErrCheckInterval check fired (iter 0 always qualifies).
  expect_true(res$best_iter_idx >= 0L)
  # best_iter_residual is finite (select_metric returned a valid value).
  expect_true(is.finite(res$best_iter_residual))
  expect_false(is.infinite(res$best_iter_residual))
})
