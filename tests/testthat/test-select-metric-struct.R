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
