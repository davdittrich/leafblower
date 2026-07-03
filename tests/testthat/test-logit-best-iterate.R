# ──────────────────────────────────────────────────────────────────────────────
# CR-C7/C8 (kxna.7/8): logit best-iterate + Armijo correctness.
#
# CR-C7: best-iterate (w_best/best_resid) was updated only at the iteration top
# from the PRE-step residual (max_b). The fresh POST-step residual (max_abs_resid)
# fed only the OK gate, so on BUDGET/STALL exits the FINAL accepted Newton step's
# improvement was discarded (logit returns w_best) and the lagged best_resid
# biased the BUDGET-vs-STALL classification. Fix: also fold the post-step residual
# into best-iterate tracking — the returned best iterate now advances one step.
#
# CR-C8: the Armijo line search halved alpha on its final failing iteration too,
# so on exhaustion the code applied the 11th, never-evaluated (2^-10-further-halved)
# step. Fix: track the best TESTED trial and, on exhaustion, apply that argmin
# step (an evaluated one), never the untested further-halved alpha.
# ──────────────────────────────────────────────────────────────────────────────

test_that("logit best-iterate captures the final post-step at BUDGET (CR-C7)", {
  set.seed(7); n <- 800L
  df <- data.frame(
    a = factor(sample(c("A", "B", "C"), n, TRUE, prob = c(.5, .3, .2))),
    b = factor(sample(c("X", "Y"),      n, TRUE, prob = c(.6, .4)))
  )
  tg <- list(a = c(A = .34, B = .33, C = .33), b = c(X = .5, Y = .5))

  # Two Newton steps (BUDGET at iter 2). Pre-fix best_error lagged one step behind
  # and reported the residual after ONE step (~5.0e-2); with the post-step update it
  # reflects BOTH steps. CR-C7b (dtk8): best_error is now measured on the FINALIZED
  # RETURNED weights (Σw=n), not the un-finalized cell iterate, so it reports the
  # residual the output actually realizes (~4.85e-3) rather than the cell-iterate
  # residual the old code returned (~6.0e-3) — returned≡reported (mxcl.5).
  r2 <- attr(suppressWarnings(harvest(df, tg, method = "logit", max_iterations = 2L,
                                      attach_weights = FALSE)), "result")
  r3 <- attr(suppressWarnings(harvest(df, tg, method = "logit", max_iterations = 3L,
                                      attach_weights = FALSE)), "result")
  expect_lt(r2$best_error, 1e-2)                 # pre-fix ~5.0e-2 → fails
  # best-iterate advanced one full step; value is the returned-weights residual (dtk8).
  expect_equal(r2$best_error, 4.852730e-03, tolerance = 1e-4)
  expect_lt(r3$best_error, r2$best_error)        # strictly better with one more step
})

test_that("logit still converges cleanly with the Armijo best-tested step (CR-C8 no-regression)", {
  # Well-conditioned feasible problem: Armijo does not exhaust, the fix is inert,
  # and logit must still reach the optimum with Σw = n.
  set.seed(11); n <- 1000L
  df <- data.frame(a = factor(sample(c("A", "B", "C"), n, TRUE)),
                   b = factor(sample(c("X", "Y"), n, TRUE)))
  tg <- list(a = c(A = 1/3, B = 1/3, C = 1/3), b = c(X = .5, Y = .5))
  w <- suppressWarnings(harvest(df, tg, method = "logit",
                                max_iterations = 100L, attach_weights = FALSE))
  r <- attr(w, "result")
  expect_identical(r$status, 0L)
  expect_lt(r$max_error, 1e-6)
  expect_equal(sum(w), n, tolerance = 1e-6)
})

test_that("logit enforces Sigma w = n and returned == reported best_error on BUDGET (CR-C7b/dtk8)", {
  # On BUDGET/STALL exit logit returned its cell iterate unscaled (Sigma w = 802 for
  # n=800) and reported the cell-iterate residual, which differs from what the
  # obs-expanded weights realize. Fix routes best_weights through finalize_weights_buf
  # (Sigma w = n) and computes best_error from the finalized RETURNED weights.
  set.seed(7L); n <- 800L
  df <- data.frame(a = factor(sample(c("x", "y", "z"), n, TRUE, prob = c(.5, .3, .2))),
                   b = factor(sample(c("s", "t"), n, TRUE, prob = c(.6, .4))))
  tg <- list(a = c(x = .34, y = .33, z = .33), b = c(s = .5, t = .5))
  o <- suppressWarnings(harvest(df, tg, method = "logit", max_iterations = 2L,
                                attach_weights = FALSE))
  r <- attr(o, "result"); w <- as.numeric(o)
  expect_equal(r$status, 4L)                       # BUDGET
  expect_equal(sum(w), n, tolerance = 1e-9)        # Sigma w = n enforced (was 802.67)

  # Recompute the max absolute-count margin residual from the RETURNED weights.
  resid <- 0
  for (v in names(tg)) {
    agg <- tapply(w, df[[v]], sum)
    for (lv in names(tg[[v]]))
      resid <- max(resid, abs(tg[[v]][[lv]] * n - agg[[lv]]) / n)
  }
  expect_equal(r$max_error, resid, tolerance = 1e-6)   # returned == reported (mxcl.5)
})
