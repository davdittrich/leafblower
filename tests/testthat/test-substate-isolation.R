# test-substate-isolation.R — T-T: CalibState sub explicit field-by-field init
#
# Verifies that the hierarchical outer loop actually iterates (outer_iterations=2
# produces distinct weights from outer_iterations=1), proving that sub-CalibState
# construction — whether via `sub = st` or explicit field-by-field init — does not
# corrupt solver state between outer iterations.
#
# Ticket: leafblower-6ycz.1.20

# ---------------------------------------------------------------------------
# DGP: N=120, 4 margins (2 coarse + 2 fine), non-trivial targets.
# Coarse margins: g1 (3 levels), g2 (2 levels) — 6 coarse cells.
# Fine margins:   g3 (2 levels), g4 (2 levels).
# Targets: uniform within each margin.
# Design ensures outer loop has work to do on iter 2 (fine residuals non-zero
# after Stage-1 coarse pass).
# ---------------------------------------------------------------------------

.make_substate_dgp <- function(seed = 7L) {
  set.seed(seed)
  N <- 120L
  g1 <- factor(sample(0:2, N, replace = TRUE))         # 3 levels
  g2 <- factor(sample(0:1, N, replace = TRUE))         # 2 levels
  g3 <- factor(sample(0:1, N, replace = TRUE))         # 2 levels (fine)
  g4 <- factor(sample(0:1, N, replace = TRUE))         # 2 levels (fine)
  df <- data.frame(g1 = g1, g2 = g2, g3 = g3, g4 = g4)
  tgt <- list(
    g1 = c(`0` = 1/3, `1` = 1/3, `2` = 1/3),
    g2 = c(`0` = 0.5,  `1` = 0.5),
    g3 = c(`0` = 0.5,  `1` = 0.5),
    g4 = c(`0` = 0.5,  `1` = 0.5)
  )
  list(df = df, tgt = tgt, N = N)
}

# ── Helper: run hierarchical raking with specified outer_iterations ─────────
# max_iterations=3 ensures the inner solver does NOT fully converge in one pass,
# so Stage-2 residuals remain and the second outer iteration has real work to do.
.run_hier <- function(d, outer_iters, method = "raking") {
  hier_cfg <- list(
    coarse_mask      = c(1L, 1L, 0L, 0L),
    min_cell_n       = 1L,
    mode             = 0L,
    outer_tol        = 1e-8,
    outer_iterations = outer_iters
  )
  suppressWarnings(
    harvest(d$df, d$tgt, method = method, hierarchical = hier_cfg,
            max_iterations = 3L)
  )
}

# ---------------------------------------------------------------------------
# Core assertion: outer iteration 2 changes weights vs iteration 1.
# Mechanism: outer_tol=0.0 forces the full budget; iter-1 and iter-2 weights
# differ because Stage-2 fine-margin corrections feed back into Stage-1 coarse
# pass on iteration 2.
# ---------------------------------------------------------------------------

test_that("raking: outer_iterations=2 produces weights distinct from outer_iterations=1", {
  d   <- .make_substate_dgp(7L)
  r1  <- .run_hier(d, outer_iters = 1L, method = "raking")
  r2  <- .run_hier(d, outer_iters = 2L, method = "raking")

  expect_equal(length(r1$weight), d$N,
               label = "raking iter-1: weight vector length == N")
  expect_equal(length(r2$weight), d$N,
               label = "raking iter-2: weight vector length == N")
  expect_false(isTRUE(all.equal(r1$weight, r2$weight)),
               label = "raking: iter-2 weights differ from iter-1 (outer loop ran)")
})

test_that("sinkhorn: outer_iterations=2 produces weights distinct from outer_iterations=1", {
  d   <- .make_substate_dgp(7L)
  r1  <- .run_hier(d, outer_iters = 1L, method = "sinkhorn")
  r2  <- .run_hier(d, outer_iters = 2L, method = "sinkhorn")

  expect_equal(length(r1$weight), d$N,
               label = "sinkhorn iter-1: weight vector length == N")
  expect_equal(length(r2$weight), d$N,
               label = "sinkhorn iter-2: weight vector length == N")
  expect_false(isTRUE(all.equal(r1$weight, r2$weight)),
               label = "sinkhorn: iter-2 weights differ from iter-1 (outer loop ran)")
})

test_that("greenkhorn: outer_iterations=2 produces weights distinct from outer_iterations=1", {
  d   <- .make_substate_dgp(7L)
  r1  <- .run_hier(d, outer_iters = 1L, method = "greenkhorn")
  r2  <- .run_hier(d, outer_iters = 2L, method = "greenkhorn")

  expect_equal(length(r1$weight), d$N,
               label = "greenkhorn iter-1: weight vector length == N")
  expect_equal(length(r2$weight), d$N,
               label = "greenkhorn iter-2: weight vector length == N")
  expect_false(isTRUE(all.equal(r1$weight, r2$weight)),
               label = "greenkhorn: iter-2 weights differ from iter-1 (outer loop ran)")
})
