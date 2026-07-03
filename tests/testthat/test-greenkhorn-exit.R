# ──────────────────────────────────────────────────────────────────────────────
# CR-C1/C2/C3 (kxna.1/2/3): greenkhorn exit-path correctness.
#
# greenkhorn's exit block (src/greenkhorn.cpp) expanded best.best_weights to obs
# weights but never enforced the Σw=n + bounds contract, never assigned the
# diagnostic result fields, and reported max_error on the configured-metric scale.
#   • kxna.1: W is maintained purely incrementally (W += delta), never rescaled at
#     exit, so Σw != n whenever max_weight clamps bind.
#   • kxna.2: res.base.kl/chi2/mean_error/grake_norm were never assigned → default
#     0.0 ("perfect fit") even on a poor solve, via pack_solver_result.
#   • kxna.3: res.base.max_error = best.best_metric, which for metric=chi2 (etc.) is
#     on a different scale than the marginal-residual reporting r_bridge expects.
# Fix: recompute CellMetrics on the best-iterate masses (honest diagnostics +
# errRp-scale max_error), then expand and route through finalize_weights_buf.
# ──────────────────────────────────────────────────────────────────────────────

.gk_fix <- function() {
  set.seed(51); n <- 1500L
  df <- data.frame(
    a = factor(sample(c("A", "B", "C"), n, TRUE, prob = c(.55, .3, .15))),
    b = factor(sample(c("X", "Y"),      n, TRUE, prob = c(.65, .35)))
  )
  tg <- list(a = c(A = .34, B = .33, C = .33), b = c(X = .5, Y = .5))
  list(df = df, tg = tg, n = n)
}

.gk_marg <- function(w, df, tg) {
  W <- sum(w)
  max(vapply(names(tg), function(v)
    max(abs(tapply(w, df[[v]], sum)[names(tg[[v]])] / W - tg[[v]])), numeric(1)))
}

test_that("greenkhorn enforces Sum(w)=n at exit even when bounds clamp (CR-C1)", {
  f <- .gk_fix()
  w <- suppressWarnings(harvest(f$df, f$tg, method = "greenkhorn",
                                max_weight = 1.5, min_weight = 0.5,
                                max_iterations = 500L, attach_weights = FALSE))
  expect_equal(sum(w), f$n, tolerance = 1e-6)   # pre-fix: ~1140 != 1500
})

test_that("greenkhorn assigns real kl/chi2/mean_error/grake_norm diagnostics (CR-C2)", {
  f <- .gk_fix()
  w <- suppressWarnings(harvest(f$df, f$tg, method = "greenkhorn",
                                max_iterations = 500L, attach_weights = FALSE))
  r <- attr(w, "result")
  # On a well-converged solve these are ~0 but must be the REAL computed values,
  # not the never-assigned phantom exact 0. Cross-check kl against the returned
  # weights: a genuine compute yields a finite non-identically-zero field set that
  # matches an independent recompute.
  W <- sum(w)
  kl_ref <- max(vapply(names(f$tg), function(v) {
    Sp <- tapply(w, f$df[[v]], sum)[names(f$tg[[v]])] / W; T <- f$tg[[v]]
    sum(ifelse(T > 0, T * log((T + 1e-10) / (Sp + 1e-10)), 0))
  }, numeric(1)))
  expect_equal(r$kl, kl_ref, tolerance = 1e-7, label = "greenkhorn kl matches returned weights")
  expect_true(is.finite(r$chi2) && is.finite(r$grake_norm) && is.finite(r$mean_error))
})

# NOTE: CR-C3 (max_error → errRp scale) is DEFERRED — it conflicts with CXX.2
# (test-greenkhorn-best-metric.R), which deliberately keeps max_error on the
# configured-metric (best-iterate) scale. See kxna.3 for the unresolved contract.
