# CR-F11 (dtkn.11): the single-pass (K>=3) and K-pass (K<3 / NA) margin_kl paths
# diverged at a targeted factor level with zero observations and a numerically-zero
# target (<= 1e-12, the only case that survives solver feasibility): single-pass
# returned NA (a tapply phantom level with W=NA defeated the Inf-guard) while K-pass
# returned Inf (droplevels + Inf-guard). Both now agree: a <= 1e-12 target on a
# zero-obs level is negligible, so margin_kl stays finite; a > 1e-12 target on a
# zero-obs level is genuine infeasibility and still returns Inf.

cqm <- leafblower:::compute_quality_metrics

test_that("single-pass and K-pass margin_kl are both finite at a zero-obs/zero-target level (CR-F11)", {
  set.seed(1); n <- 60
  w <- rep(1, n)

  # Single-pass fixture: K=3 factors, no NA; margin 'a' level 'z' has 0 obs, target 1e-12.
  df3 <- data.frame(
    a = factor(sample(c("x", "y"), n, TRUE), levels = c("x", "y", "z")),
    b = factor(sample(c("p", "q"), n, TRUE)),
    c = factor(sample(c("m", "n"), n, TRUE))
  )
  tg3 <- list(a = c(x = 0.5, y = 0.5 - 1e-12, z = 1e-12),
              b = c(p = 0.5, q = 0.5),
              c = c(m = 0.5, n = 0.5))
  mk_sp <- cqm(w, tg3, df3)$margin_kl

  # K-pass fixture: K=2 factors; margin 'a' level 'z' has 0 obs, target 1e-12.
  df2 <- data.frame(
    a = factor(sample(c("x", "y"), n, TRUE), levels = c("x", "y", "z")),
    b = factor(sample(c("p", "q"), n, TRUE))
  )
  tg2 <- list(a = c(x = 0.5, y = 0.5 - 1e-12, z = 1e-12),
              b = c(p = 0.5, q = 0.5))
  mk_kp <- cqm(w, tg2, df2)$margin_kl

  expect_true(is.finite(mk_sp))  # was NA pre-fix
  expect_true(is.finite(mk_kp))  # was Inf pre-fix
})

test_that("a > 1e-12 target on a zero-obs level is still Inf (genuine infeasibility, CR-F11)", {
  set.seed(2); n <- 60
  w <- rep(1, n)
  # K-pass: 'z' has 0 obs but a non-negligible target 0.2 → unmatched margin → Inf.
  df2 <- data.frame(
    a = factor(sample(c("x", "y"), n, TRUE), levels = c("x", "y", "z")),
    b = factor(sample(c("p", "q"), n, TRUE))
  )
  tg2 <- list(a = c(x = 0.4, y = 0.4, z = 0.2), b = c(p = 0.5, q = 0.5))
  expect_true(is.infinite(cqm(w, tg2, df2)$margin_kl))
})
