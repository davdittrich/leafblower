# ──────────────────────────────────────────────────────────────────────────────
# CR-F4/F5 (dtkn.4/5): diagnostic-path guards.
#
# CR-F4: harvest.R's BUDGET warning used base `%||%` (R>=4.4 only, no version
# floor declared) on e_final, which is never NULL there (REALSXP scalar) — a dead
# guard that ALSO crashed the diagnostic with "could not find function '%||%'" on
# R<4.4. The `%||% NaN` was removed. Regression guard: the BUDGET warning fires.
#
# CR-F5: get_current_miss() did `col <- data[[v]]` with no NULL check, so an absent
# variable silently produced W=0/prop=0 and returned max(tgt) as a fabricated
# "miss". Now stop()s, matching diagnose_weights().
# ──────────────────────────────────────────────────────────────────────────────

test_that("get_current_miss stops on an absent variable, not fabricate a miss (CR-F5)", {
  set.seed(1); n <- 200L
  data <- data.frame(a = factor(sample(c("A", "B"), n, TRUE)))
  target <- list(absent_var = c(A = 0.5, B = 0.5))
  weights <- rep(1, n)
  expect_error(get_current_miss(data, target, weights),
               "not found in data")
})

test_that("get_current_miss works for present variables (CR-F5 no-regression)", {
  set.seed(1); n <- 200L
  data <- data.frame(a = factor(sample(c("A", "B"), n, TRUE)))
  target <- list(a = c(A = 0.5, B = 0.5))
  m <- get_current_miss(data, target, rep(1, n))
  expect_true(is.finite(m[["a"]]) && m[["a"]] >= 0)
})

test_that("BUDGET-exhausted warning path executes without a %||% error (CR-F4)", {
  set.seed(7); n <- 800L
  df <- data.frame(
    a = factor(sample(c("A", "B", "C"), n, TRUE, prob = c(.5, .3, .2))),
    b = factor(sample(c("X", "Y"),      n, TRUE, prob = c(.6, .4)))
  )
  tg <- list(a = c(A = .34, B = .33, C = .33), b = c(X = .5, Y = .5))
  # max_iterations=2 forces a BUDGET(4) exit → the warning path at harvest.R:700-708.
  expect_warning(
    harvest(df, tg, method = "logit", max_iterations = 2L, attach_weights = FALSE),
    "budget exhausted"
  )
})
