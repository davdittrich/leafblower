test_that("all solvers produce same output after select_metric struct migration", {
  set.seed(99)
  n   <- 60L
  df  <- data.frame(
    x = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
    y = factor(sample(c("p", "q"),       n, replace = TRUE))
  )
  pop <- list(x = c(a = 1/3, b = 1/3, c = 1/3), y = c(p = 0.5, q = 0.5))

  for (m in c("oris", "chebyshev", "raking", "sinkhorn")) {
    w <- suppressWarnings(
      leafblower::harvest(df, pop, method = m, max_iterations = 500L,
                          convergence = list(absolute = 1e-3))
    )
    r <- attr(w, "result")
    expect_lte(r$max_error, 1e-3, label = paste("max_error for", m))
  }
})
# ──────────────────────────────────────────────────────────────────────────────
# CR-A3 (mxcl.3): metric="marginal_kl" must drive REAL convergence, not a phantom
# zero. The convenience select_metric(CalibMetric, CellMetrics) overload omitted
# marginal_kl (defaulted to 0.0), so check_convergence instant-triggered — greenkhorn
# with metric="marginal_kl" stopped at the first check interval (iter 10, achieved
# ~1.9e-4) instead of converging like metric="kl" (iter 30, ~6.6e-11). Fix: track
# a true Σ_k marginal_kl in CellMetrics and pass it through select_metric.
# ──────────────────────────────────────────────────────────────────────────────

test_that("metric=marginal_kl converges properly, not spurious phantom-zero (CR-A3)", {
  set.seed(707)
  n <- 2000L
  a <- factor(sample(c("A", "B", "C"), n, replace = TRUE, prob = c(0.6, 0.25, 0.15)))
  b <- factor(sample(c("X", "Y"),      n, replace = TRUE, prob = c(0.7, 0.3)))
  data <- data.frame(a = a, b = b)
  target <- list(a = c(A = 0.33, B = 0.34, C = 0.33), b = c(X = 0.5, Y = 0.5))

  achieved <- function(w) {
    Wt <- sum(w)
    max(vapply(names(target), function(v)
      max(abs(tapply(w, data[[v]], sum)[names(target[[v]])] / Wt - target[[v]])),
      numeric(1)))
  }

  # All solvers affected by the phantom-zero (raking hand-builds its CellMetrics,
  # so it needed a per-solver fix too — raking was the worst: instant iter-1 OK
  # at achieved ~4.9e-3; greenkhorn/sinkhorn/chebyshev consume compute_cell_metrics).
  for (m in c("raking", "greenkhorn", "sinkhorn", "chebyshev")) {
    w_mkl <- suppressWarnings(harvest(data, target, method = m,
                                      convergence = list(metric = "marginal_kl"),
                                      max_iterations = 500L, attach_weights = FALSE, verbose = 0))
    w_kl <- suppressWarnings(harvest(data, target, method = m,
                                     convergence = list(metric = "kl"),
                                     max_iterations = 500L, attach_weights = FALSE, verbose = 0))
    # marginal_kl must actually calibrate (not phantom-stop) and match the honest
    # kl metric's fixed point.
    expect_lt(achieved(w_mkl), 1e-8, label = paste("marginal_kl achieved for", m))
    expect_equal(achieved(w_mkl), achieved(w_kl), tolerance = 1e-9,
                 label = paste("marginal_kl vs kl for", m))
  }
})
