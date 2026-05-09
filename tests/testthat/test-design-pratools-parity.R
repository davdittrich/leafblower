test_that("design_effect 4-arg matches PracTools::deffH under uniform weights (rtol=1e-10)", {
  skip_if_not_installed("PracTools")
  set.seed(2024L); n <- 200L
  data <- data.frame(
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  target <- list(region = c(N = 0.25, S = 0.25, E = 0.25, W = 0.25))
  w <- rep(1.0, n)  # uniform weights -> var(w)=0 -> PracTools Eq3.4 second-component vanishes
  y <- 10 + 3 * (data$region == "N") - 2 * (data$region == "S") + rnorm(n)

  # leafblower output (Eq 3.5 via C++ core)
  d_lbw <- design_effect(w, outcome = y, data = data, target = target)

  # PracTools direct oracle (Eq 3.4; under uniform w, reduces to Eq 3.5)
  # X uses target level order, dropping first level (no intercept) — same as C++ build
  X <- model.matrix(~ region, data = data)[, -1L, drop = FALSE]
  d_pratools <- PracTools::deffH(w = w, y = y, x = X)

  expect_equal(d_lbw, d_pratools, tolerance = 1e-10,
               label = sprintf("d_lbw=%.15f vs PracTools::deffH=%.15f (uniform-w, rtol=1e-10)",
                               d_lbw, d_pratools))
})

# Survey fallback oracle: hand-roll Eq 3.5 matching C++ arithmetic exactly.
# Levels taken from target$region names (N S E W); code 0 = N dropped → columns S, E, W.
# Tolerance 1e-6: floating-point accumulation across n=200 obs differs by ~4e-8 (relative)
# between R and C++ due to operation order; rtol=1e-10 is not achievable at this n.
test_that("design_effect 4-arg matches hand-rolled Eq 3.5 oracle (rtol=1e-6, survey fallback)", {
  skip_if_not_installed("survey")
  set.seed(2024L); n <- 200L
  data2 <- data.frame(
    region = sample(c("N", "S", "E", "W"), n, replace = TRUE),
    stringsAsFactors = FALSE
  )
  target2 <- list(region = c(N = 0.25, S = 0.25, E = 0.25, W = 0.25))
  w2 <- rep(1.0, n)
  y2 <- 10 + 3 * (data2$region == "N") - 2 * (data2$region == "S") + rnorm(n)
  d_lbw2 <- design_effect(w2, outcome = y2, data = data2, target = target2)

  # Reproduce C++ X matrix: levels from target names, 0-based codes, drop level 0.
  lvls <- names(target2$region)           # N S E W
  codes <- as.integer(factor(data2$region, levels = lvls)) - 1L  # 0=N,1=S,2=E,3=W
  p <- length(lvls) - 1L                 # 3 columns: S E W
  X2 <- matrix(0.0, nrow = n, ncol = p)
  for (i in seq_len(n)) { c_i <- codes[i]; if (c_i >= 1L) X2[i, c_i] <- 1.0 }

  # WLS beta via Cholesky (matches C++ dpotrf path up to float order)
  XtWX <- t(X2) %*% (w2 * X2)
  XtWy <- as.vector(t(X2) %*% (w2 * y2))
  L <- chol(XtWX)
  beta_hat <- backsolve(L, forwardsolve(t(L), XtWy))
  u2 <- y2 - as.vector(X2 %*% beta_hat)

  wtd_var <- function(z, w) {
    w_bar <- sum(w * z) / sum(w)
    sum(w * (z - w_bar)^2) / sum(w)
  }
  deff_K2 <- length(w2) * sum(w2^2) / sum(w2)^2
  d_analytic2 <- deff_K2 * wtd_var(u2, w2) / wtd_var(y2, w2)

  expect_equal(d_lbw2, d_analytic2, tolerance = 1e-6,
               label = sprintf("d_lbw=%.15f vs analytic=%.15f (rel diff=%.2e)",
                               d_lbw2, d_analytic2,
                               abs(d_lbw2 - d_analytic2) / d_analytic2))
})
