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
  # The improvement/pct rule fires on a metric PLATEAU; on a bounds-blocked problem
  # (tight max/min_weight + targets far from the empirical margins) the metric
  # freezes from iter 1 because weights are pinned at their bounds, so the OLD code
  # marked RK_OK despite grossly unmet margins. Fix (raking.cpp kRakingFeasTol=1%):
  # RK_OK requires errRp <= 1% feasibility; a plateau above that -> STALL(5).
  set.seed(88); n <- 1000L
  df <- data.frame(
    a = factor(sample(c("p", "q", "r"), n, TRUE, prob = c(.7, .2, .1))),
    b = factor(sample(c("s", "t"),      n, TRUE, prob = c(.6, .4)))
  )
  tg <- list(a = c(p = .2, q = .3, r = .5),   # far from empirical (.7,.2,.1)
             b = c(s = .3, t = .7))           # far from empirical (.6,.4)
  w <- suppressWarnings(harvest(df, tg, method = "raking",
                                max_weight = 1.3, min_weight = 0.7,
                                max_iterations = 500L, attach_weights = FALSE))
  r <- attr(w, "result")
  expect_gt(r$max_error, 0.1)     # genuinely bounds-blocked: margins grossly unmet
  expect_false(r$status == 0L)    # must NOT be RK_OK
  expect_equal(r$status, 5L)      # STALL: at constrained optimum, weights valid
})

test_that("raking feasible convergence still reports RK_OK (jy0m gate not over-strict)", {
  # Guard the fix against over-strictness: a feasible problem whose margins are met
  # (errRp well inside the 1pp bar) must still return RK_OK, not STALL.
  set.seed(1); n <- 1000L
  df <- data.frame(a = factor(sample(c("p", "q", "r"), n, TRUE)),
                   b = factor(sample(c("s", "t"), n, TRUE)))
  tg <- list(a = c(p = 1/3, q = 1/3, r = 1/3), b = c(s = .5, t = .5))
  w <- suppressWarnings(harvest(df, tg, method = "raking",
                                max_weight = 5, min_weight = 0.2,
                                max_iterations = 500L, attach_weights = FALSE))
  r <- attr(w, "result")
  expect_lt(r$max_error, 1e-2)    # feasibly met margins (inside the bar)
  expect_equal(r$status, 0L)      # RK_OK preserved
})

test_that("raking SRAA (accelerate=TRUE) bounds-blocked plateau exits STALL(5) (jy0m.2)", {
  # The SRAA path shares the scale-blind improvement-rule plateau bug. The same
  # feasibility gate applies there; a bounds-blocked accelerated run must NOT report
  # RK_OK. jy0m.2 added the SRAA STALL emitter, so a bounds-blocked plateau now exits
  # STALL(5) (constrained optimum) instead of running to f_evals exhaustion → BUDGET(4).
  set.seed(88); n <- 1000L
  df <- data.frame(a = factor(sample(c("p", "q", "r"), n, TRUE, prob = c(.7, .2, .1))),
                   b = factor(sample(c("s", "t"), n, TRUE, prob = c(.6, .4))))
  tg <- list(a = c(p = .2, q = .3, r = .5), b = c(s = .3, t = .7))
  w <- suppressWarnings(harvest(df, tg, method = "raking", accelerate = TRUE,
                                max_weight = 1.3, min_weight = 0.7,
                                max_iterations = 500L, attach_weights = FALSE))
  r <- attr(w, "result")
  expect_gt(r$max_error, 0.1)
  expect_equal(r$status, 5L)      # STALL(5): constrained optimum, not false RK_OK, not BUDGET
})

test_that("box-infeasibility overrides even an explicit user tolerance → STALL(5) (jy0m.1)", {
  # An EXPLICIT loose user tolerance (convergence=list(absolute=0.6)) accepts errRp up
  # to 0.6, and the seed-88 fixture plateaus at errRp≈0.5. But this fixture is PROVABLY
  # box-infeasible (a=p: target .2 vs empirical .7 needs ratio 0.29 < min_weight 0.7),
  # so RK_OK would be a false success claim on an unreachable margin. The box
  # certificate denies RK_OK regardless of the user's tolerance → STALL(5); the
  # best-effort weights are still returned. (Feasible problems above 1pp still honor
  # the user tol — see the min_weight=0.2 over-strictness test above, which is
  # box-feasible.)
  set.seed(88); n <- 1000L
  df <- data.frame(a = factor(sample(c("p", "q", "r"), n, TRUE, prob = c(.7, .2, .1))),
                   b = factor(sample(c("s", "t"), n, TRUE, prob = c(.6, .4))))
  tg <- list(a = c(p = .2, q = .3, r = .5), b = c(s = .3, t = .7))
  w <- suppressWarnings(harvest(df, tg, method = "raking",
                                max_weight = 1.3, min_weight = 0.7, max_iterations = 500L,
                                attach_weights = FALSE,
                                convergence = list(absolute = 0.6, metric = "max_err")))
  r <- attr(w, "result")
  expect_equal(r$status, 5L)      # box-infeasible ⇒ STALL(5), even under a loose user tol
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

test_that("raking SRAA metric='l1_weight' + abs tol not bypassed by phantom l1 (jy0m.3)", {
  # compute_cell_metrics never writes CellMetrics.l1 (stays 0). Before jy0m.3, the SRAA
  # gate's user_abs_met = select_metric(l1)=0 < absolute_tol was spuriously true from
  # iter 1, marking a bounds-blocked run RK_OK regardless of margins. jy0m.3 populates
  # last_F_metrics.l1 with the honest obs-level Σ|ΔX|/n BEFORE the select_metric read,
  # so the feasibility gate (jy0m.1 box-cert + errRp) is no longer bypassed and the
  # box-infeasible seed-88 fixture correctly exits STALL(5), not a false RK_OK(0).
  set.seed(88); n <- 1000L
  df <- data.frame(a = factor(sample(c("p", "q", "r"), n, TRUE, prob = c(.7, .2, .1))),
                   b = factor(sample(c("s", "t"), n, TRUE, prob = c(.6, .4))))
  tg <- list(a = c(p = .2, q = .3, r = .5), b = c(s = .3, t = .7))
  w <- suppressWarnings(harvest(df, tg, method = "raking", accelerate = TRUE,
                                max_weight = 1.3, min_weight = 0.7, max_iterations = 500L,
                                attach_weights = FALSE,
                                convergence = list(absolute = 1e-6, metric = "l1_weight")))
  r <- attr(w, "result")
  expect_false(r$status == 0L)    # phantom-l1 no longer yields a false RK_OK
  expect_equal(r$status, 5L)      # STALL(5) via box-cert + honest l1
})
