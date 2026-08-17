test_that("design_effect 4-arg matches PracTools::deffK/glm Eq-3.5 oracle (rtol=1e-6)", {
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

  # X uses target level order, dropping first level (no intercept) — same as C++ build
  X <- model.matrix(~ region, data = data)[, -1L, drop = FALSE]

  # Upstream-singularity tripwire: PracTools::deffH's Eq 3.4 correction term forms
  # rho.u2w and rho.uw as cov/sqrt(vu2*vw); under uniform w, vw = wtdvar(w, w) = 0, so
  # both ratios are 0/0 = NaN and the term evaluates 0*NaN = NaN instead of vanishing
  # (CRAN PracTools 1.7.5, R/deffH.R:22-25 — a version-independent removable
  # singularity, not something a version pin can fix). Only compare when finite, so the
  # direct third-party call stays in the suite and self-heals if upstream ever removes
  # the singularity, without ever asserting on a NaN.
  d_full <- PracTools::deffH(w = w, y = y, x = X)
  if (is.finite(d_full)) {
    expect_equal(d_lbw, d_full, tolerance = 1e-6)
  }

  # Eq-3.5 oracle built from primitives that ARE defined at var(w) = 0: PracTools::deffK
  # for the (third-party) Kish factor, and a WLS fit via glm() — a QR path fully
  # independent of the C++ core's Cholesky path — for u = coef[1] + residuals, exactly
  # as PracTools forms u internally.
  reg <- glm(y ~ X, weights = w)
  u <- unname(reg$coefficients[1L]) + reg$residuals
  # Local weighted variance. PracTools::wtdvar carries an n/(n-1) factor that cancels in
  # this ratio, so this local helper and PracTools' own helper give the same oracle.
  wv <- function(z, wt) {
    m <- sum(wt * z) / sum(wt)
    sum(wt * (z - m)^2) / sum(wt)
  }
  d_oracle <- PracTools::deffK(w) * wv(u, w) / wv(y, w)

  expect_equal(d_lbw, d_oracle, tolerance = 1e-6,
               label = sprintf("d_lbw=%.15f vs oracle=%.15f (PracTools::deffK * glm-WLS Eq 3.5)",
                               d_lbw, d_oracle))
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
