context("oris faithful — non-uniform design weights within cell")

test_that("marginals hit targets when d[i] varies within cell", {
  set.seed(42)
  n <- 1000
  # Only 2 cells total: (a=0, b=0) and (a=1, b=0). Force wide d[i] variation within each.
  # Columns MUST be factor or character (r_bridge.cpp:145 errors on numeric).
  df <- data.frame(
    a = factor(rep(0:1, each = n/2)),
    b = factor(rep(0L, n))
  )
  # d[i] varies 1 to 10 within each cell
  d <- rep(c(1, 10), length.out = n)
  tgt <- list(
    a = c(`0` = 0.3, `1` = 0.7),
    b = c(`0` = 1.0)
  )
  # max_weight=2.0 is tight: within cell a=1 (target 0.7), d=10 obs would land
  # at w=2.55 without clamp — this is the intentional design post-fix. Pre-fix
  # clamp pinned them to 2.0 and distorted the a=1 marginal by ~7.7%.
  res <- harvest(df, tgt, method = "oris",
                 start_weights = d, max_weight = 2.0, min_weight = 0,
                 convergence = list(absolute = 1e-6))
  w <- res$weights
  # PRIMARY assertion: marginals hit targets (the bug breaks this).
  m_a0 <- sum(w[df$a == 0]) / sum(w)
  m_a1 <- sum(w[df$a == 1]) / sum(w)
  expect_lt(abs(m_a0 - 0.3), 1e-4)
  expect_lt(abs(m_a1 - 0.7), 1e-4)
  # CELL-AGGREGATE assertion: per-cell mean weight respects max_weight. This is
  # the bound the faithful ORIS guarantees (X[c] ≤ U_cell[c] = max_weight * n_per_cell[c]).
  # Per-OBS bounds intentionally can leak by d_max/d_mean factor on non-uniform d;
  # this trade-off is documented in src/oris.cpp expansion comment.
  cell_mean_a0 <- mean(w[df$a == 0])
  cell_mean_a1 <- mean(w[df$a == 1])
  expect_lte(cell_mean_a0, 2.0 + 1e-8)
  expect_lte(cell_mean_a1, 2.0 + 1e-8)
})
