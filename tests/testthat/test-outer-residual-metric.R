# test-outer-residual-metric.R
# Spec §8 coarse-margin L1 metric: reported hier_outer_residual_final ==
# spec L1 computed in R from returned weights against coarse targets only.
# Tests the BUDGET-exit recompute fix (T-O critical fix).
#
# Formula: (Σ_k Σ_j |Σ_i w_i · 1[X_ki=j] − target_kj·N|) / N
# over coarse margins only.

library(leafblower)

# ── helper: compute spec L1 in R ─────────────────────────────────────────────
spec_l1 <- function(weights, df, targets, coarse_mask) {
  N     <- length(weights)
  K     <- ncol(df)
  total <- 0.0
  for (k in seq_len(K)) {
    if (coarse_mask[k] != 1L) next
    tgt  <- targets[[k]]
    cats <- names(tgt)
    col  <- as.character(df[[k]])
    for (j in seq_along(cats)) {
      achieved <- sum(weights[col == cats[j]])
      total    <- total + abs(achieved - tgt[j] * N)
    }
  }
  total / N
}

# ── DGP: small orthogonal, budget-forcing ────────────────────────────────────
# K=6 binary margins; 3 coarse, 3 fine; N=150; outer_iterations deliberately
# tiny (2) to force BUDGET-exit on seeds that don't converge immediately.
make_budget_dgp <- function(seed, N = 150L) {
  set.seed(seed)
  g1 <- as.factor(rbinom(N, 1L, 0.3))
  g2 <- as.factor(rbinom(N, 1L, 0.4))
  g3 <- as.factor(rbinom(N, 1L, 0.5))
  # Fine margins correlated with coarse so orthogonality holds
  f1 <- as.factor(ifelse(as.integer(g1) == 1L,
                         rbinom(N, 1L, 0.1), rbinom(N, 1L, 0.9)))
  f2 <- as.factor(ifelse(as.integer(g2) == 1L,
                         rbinom(N, 1L, 0.1), rbinom(N, 1L, 0.9)))
  f3 <- as.factor(ifelse(as.integer(g3) == 1L,
                         rbinom(N, 1L, 0.1), rbinom(N, 1L, 0.9)))
  df      <- data.frame(g1, g2, g3, f1, f2, f3)
  targets <- lapply(df, function(col) {
    lvls <- levels(col)
    setNames(rep(1.0 / length(lvls), length(lvls)), lvls)
  })
  coarse_mask <- c(1L, 1L, 1L, 0L, 0L, 0L)
  list(df = df, targets = targets, coarse_mask = coarse_mask)
}

# ── core assertion ────────────────────────────────────────────────────────────
check_residual_matches <- function(dgp, solver, outer_iterations = 2L,
                                   tol = 1e-10, label = "") {
  res <- harvest(
    data         = dgp$df,
    target       = dgp$targets,
    method       = solver,
    hierarchical = list(
      coarse_mask      = dgp$coarse_mask,
      min_cell_n       = 1L,
      mode             = 0L,        # Strategy A
      outer_tol        = 1e-12,      # tight enough that 2 iters won't converge → BUDGET-exit
      outer_iterations = outer_iterations
    )
  )
  w        <- res$weights
  reported <- attr(res, "result")$outer_residual_final
  expected <- spec_l1(w, dgp$df, dgp$targets, dgp$coarse_mask)
  diff     <- abs(reported - expected)
  expect_lt(diff, tol,
    label = paste0(label, "solver=", solver,
                   " reported=", signif(reported, 6),
                   " expected=", signif(expected, 6)))
}

# ── run 3 seeds × 3 solvers ──────────────────────────────────────────────────
seeds   <- c(1L, 42L, 123L)
solvers <- c("raking", "sinkhorn", "greenkhorn")

for (seed in seeds) {
  dgp <- make_budget_dgp(seed)
  for (solver in solvers) {
    check_residual_matches(dgp, solver,
                           label = paste0("seed=", seed, " "))
  }
}
