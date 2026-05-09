# Tests for R/design_effect.R: thin wrapper over .Call(C_rk_design_effect, ...).

test_that("design_effect 1-arg matches Kish (1965) formula", {
  w <- c(1, 2, 3, 4)
  expected <- length(w) * sum(w^2) / sum(w)^2
  expect_equal(design_effect(w), expected, tolerance = 1e-12)
})

test_that("effective_sample_size = n / design_effect (1-arg path)", {
  w <- c(1, 2, 3, 4)
  expect_equal(effective_sample_size(w), length(w) / design_effect(w), tolerance = 1e-12)
})

test_that("design_effect 4-arg empty target -> deff_K (degenerate)", {
  w <- c(1, 2, 3, 4)
  y <- c(10, 20, 30, 40)
  data <- data.frame(g = c("a", "b", "a", "b"), stringsAsFactors = FALSE)
  expect_equal(design_effect(w, outcome = y, data = data, target = list()),
               design_effect(w), tolerance = 1e-12)
})

test_that("design_effect 4-arg 3-level perfect-fit -> deff_H = 1/3", {
  # K=1, 3 levels (A/B/C, 20 obs each), uniform w=1, y = group mean (perfect fit).
  # X = [ind_B, ind_C] (no intercept, drop-first encoding).
  # beta_hat = [20, 30]; residuals: A-group u=10, B/C u=0.
  # var_u = mean((u - u_bar)^2) = 200/9; var_y = mean((y - y_bar)^2) = 200/3.
  # deff_H = deff_K * var_u / var_y = 1 * (1/3) = 1/3.
  n <- 60L
  g <- rep(c("A", "B", "C"), each = 20L)
  y <- c(rep(10.0, 20L), rep(20.0, 20L), rep(30.0, 20L))
  data   <- data.frame(g = g, stringsAsFactors = FALSE)
  target <- list(g = c(A = 1/3, B = 1/3, C = 1/3))
  d <- design_effect(rep(1.0, n), outcome = y, data = data, target = target)
  expect_equal(d, 1/3, tolerance = 1e-6)
})
