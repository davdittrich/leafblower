# ──────────────────────────────────────────────────────────────────────────────
# CR-C5/C6 (kxna.5/6): raking exit-path correctness.
#
# CR-C5: removed the `metric==KL && wkl_flat <= tol_abs` RK_OK shortcut in
# raking.cpp. wkl_flat = Σ X·log(X/X_init)/n measures distance FROM the design
# weights, not the achieved-margin KL, so it was an incorrect OK criterion. In
# practice it is near-unreachable (check_convergence at :495 fires first), so its
# removal is safe cleanup with no observable behaviour change — hence no isolating
# regression test. The OBSERVABLE false-OK on bounds-blocked infeasible plateaus
# (status=0 at kl=0.637) is a DISTINCT improvement-rule scale-blindness bug, filed
# as leafblower-jy0m (CR-C5b); do not conflate it with the shortcut here.
#
# CR-C6: raking's SRAA (accelerate=TRUE) loop assigned only max_error/iterations;
# mean_error/kl/grake_norm/l1_weight_change were written only on the flat path, so
# accelerate=TRUE returned them as zeros ("perfect fit"), misrepresenting quality.
# ──────────────────────────────────────────────────────────────────────────────

test_that("raking bounds-blocked infeasible plateau must not report RK_OK (CR-C5b/jy0m)", {
  skip("tracked in leafblower-jy0m: improvement-rule feasibility gate, not the (removed) wkl shortcut")
})

test_that("raking SRAA (accelerate=TRUE) returns real diagnostics, not zeros (CR-C6)", {
  set.seed(88); n <- 1000L
  df <- data.frame(a = factor(sample(c("A", "B", "C"), n, TRUE)),
                   b = factor(sample(c("X", "Y"), n, TRUE)))
  tg <- list(a = c(A = 1/3, B = 1/3, C = 1/3), b = c(X = .5, Y = .5))
  w <- suppressWarnings(harvest(df, tg, method = "raking", accelerate = TRUE,
                                max_iterations = 500L, attach_weights = FALSE))
  r <- attr(w, "result")
  W <- sum(w)
  kl_ref <- max(vapply(names(tg), function(v) {
    Sp <- tapply(w, df[[v]], sum)[names(tg[[v]])] / W; T <- tg[[v]]
    sum(ifelse(T > 0, T * log((T + 1e-10) / (Sp + 1e-10)), 0))
  }, numeric(1)))
  expect_equal(r$kl, kl_ref, tolerance = 1e-7)
  expect_true(is.finite(r$mean_error) && is.finite(r$grake_norm) &&
              is.finite(r$l1_weight_change))
})
