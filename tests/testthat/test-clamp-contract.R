library(leafblower)

# CR-D11 (j7x8.11): unified exit contract across all cell-table solvers.
# On a bound-active CELL-mode problem the per-obs clamp used to distort
# marginals and break Sw=n for raking/greg/sinkhorn (measured 13pp drift).
# All solvers now share the canonical no-clamp cell contract: expand without
# per-obs clamp, then finalize_weights (Sw=n; cell mode counts violations only).
make_bound_active <- function() {
  set.seed(20260703)
  n <- 4000L
  df <- data.frame(
    x = factor(sample(c("a","b"), n, replace = TRUE)),
    y = factor(sample(c("p","q"), n, replace = TRUE))
  )
  dw <- rep(1.0, n)
  dw[sample(seq_len(n), 200L)] <- 8.0   # within-cell skew -> exit clamp would bind
  target <- list(x = c(a = 0.70, b = 0.30), y = c(p = 0.35, q = 0.65))
  list(df = df, dw = dw, target = target, n = n)
}

test_that("all cell-table solvers preserve Sw=n and marginals in cell mode (CR-D11)", {
  f <- make_bound_active()
  for (m in c("raking", "greg", "sinkhorn", "oris", "chebyshev")) {
    w <- as.numeric(harvest(f$df, f$target, method = m, design_weights = f$dw,
                            max_weight = 3, min_weight = 0.2, bounds_mode = "cell",
                            max_iterations = 2000L, attach_weights = FALSE))
    # (1) Sum-to-n preserved exactly (the guarantee the clamp used to break).
    expect_equal(sum(w), f$n, tolerance = 1e-6,
                 info = sprintf("%s: Sw=%.6f", m, sum(w)))
    # (2) Marginals met to solver tolerance (clamp used to drift 13pp).
    dev_x <- max(abs(tapply(w, f$df$x, sum) / f$n - f$target$x[levels(f$df$x)]))
    dev_y <- max(abs(tapply(w, f$df$y, sum) / f$n - f$target$y[levels(f$df$y)]))
    # 5e-3 threshold: 25x below the 13pp (0.13) regression signal, well clear of
    # the default convergence tol (1e-3) so a tol-boundary stop cannot flake it.
    expect_lt(max(dev_x, dev_y), 5e-3)
  }
})
