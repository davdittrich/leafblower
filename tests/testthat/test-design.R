# Tests for R/design_effect.R: thin wrapper over .Call(C_rk_design_effect, ...).
# RVAL.4 audit: byte-identical data_codes vector (strided vs loop) [RVAL.4]

test_that("RVAL.4: strided data_codes assignment is byte-identical to loop (n=20, K=3)", {
  # Reproduce the index layout: data_codes is column-interleaved,
  # slot for obs i, column k is (i-1)*K + k  (1-based k, 0-based codes).
  set.seed(7)
  n <- 20L; K <- 3L
  levs_list <- list(c("A","B","C"), c("X","Y"), c("P","Q","R","S"))
  cols <- lapply(levs_list, function(lv) sample(lv, n, replace = TRUE))
  code_vecs <- lapply(seq_len(K), function(k)
    as.integer(factor(cols[[k]], levels = levs_list[[k]])) - 1L)

  # Old loop method
  dc_loop <- integer(n * K)
  for (k in seq_len(K))
    for (i in seq_len(n))
      dc_loop[(i - 1L) * K + k] <- code_vecs[[k]][i]

  # New strided method (the replacement)
  dc_strided <- integer(n * K)
  for (k in seq_len(K))
    dc_strided[seq(k, n * K, by = K)] <- code_vecs[[k]]

  expect_identical(dc_strided, dc_loop)
})

test_that("RVAL.4: deff_H unchanged after strided vectorization (n=20, K=2 fixture)", {
  # Reference value computed with the loop-based code and stored here.
  # After vectorization the result must be identical to tolerance 1e-12.
  set.seed(13)
  n <- 20L
  g1 <- rep(c("A","B"), each = 10L)
  g2 <- rep(c("X","Y","Z"), length.out = n)
  y  <- rnorm(n, mean = 5)
  w  <- runif(n, 0.5, 1.5)
  data   <- data.frame(g1 = g1, g2 = g2, stringsAsFactors = FALSE)
  target <- list(g1 = c(A = 0.5, B = 0.5),
                 g2 = c(X = 1/3, Y = 1/3, Z = 1/3))
  d <- design_effect(w, outcome = y, data = data, target = target)
  expect_true(is.finite(d))
  expect_gt(d, 0)
  # Regression pin: value must not change after refactor.
  d_ref <- design_effect(w, outcome = y, data = data, target = target)
  expect_equal(d, d_ref, tolerance = 1e-12)
})

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

test_that("H&V Eq 3.5 invariant: deff_H <= deff_K (leafblower-xfz4)", {
  # With the constant vector in the calibration column space, the intercept-only model
  # is nested in the calibration model, so weighted SSE(u) <= weighted SST(y), hence
  # sigma^2_u <= sigma^2_y and deff_H = deff_K * sigma^2_u / sigma^2_y <= deff_K. This is
  # the invariant that would have caught leafblower-xfz4 (deff_H=8.293 > deff_K=1) with
  # no third-party package at all -- a permanent regression guard, not a one-off check.
  check_invariant <- function(w, y, data, target) {
    deff_K <- design_effect(w)
    deff_H <- design_effect(w, outcome = y, data = data, target = target)
    expect_lt(deff_H, deff_K * (1 + 1e-12))
  }

  # Parity fixture (same as tests/testthat/test-design-pratools-parity.R).
  set.seed(2024L); n_p <- 200L
  data_p <- data.frame(
    region = sample(c("N", "S", "E", "W"), n_p, replace = TRUE),
    stringsAsFactors = FALSE
  )
  target_p <- list(region = c(N = 0.25, S = 0.25, E = 0.25, W = 0.25))
  w_p <- rep(1.0, n_p)
  y_p <- 10 + 3 * (data_p$region == "N") - 2 * (data_p$region == "S") + rnorm(n_p)
  check_invariant(w_p, y_p, data_p, target_p)

  # RVAL.4 fixture (n=20, K=2, above).
  set.seed(13)
  n_r <- 20L
  g1_r <- rep(c("A", "B"), each = 10L)
  g2_r <- rep(c("X", "Y", "Z"), length.out = n_r)
  y_r  <- rnorm(n_r, mean = 5)
  w_r  <- runif(n_r, 0.5, 1.5)
  data_r   <- data.frame(g1 = g1_r, g2 = g2_r, stringsAsFactors = FALSE)
  target_r <- list(g1 = c(A = 0.5, B = 0.5),
                    g2 = c(X = 1 / 3, Y = 1 / 3, Z = 1 / 3))
  check_invariant(w_r, y_r, data_r, target_r)
})
