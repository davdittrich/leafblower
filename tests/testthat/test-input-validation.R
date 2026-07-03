# ──────────────────────────────────────────────────────────────────────────────
# CR-D2 (j7x8.2): reject non-finite bounds and invalid iteration budgets.
#
# NaN bounds slip past the `min_weight >= max_weight` check (all NaN comparisons
# are false) and NA/INT_MIN iteration budgets passed silently, both propagating to
# garbage weights returned as RK_OK. harvest.R now guards these before the C call,
# and validation.cpp / calib_validate_preentry guard both bridges' C entry.
# +Inf max_weight remains legal (unbounded).
# ──────────────────────────────────────────────────────────────────────────────

.iv_fixture <- function() {
  set.seed(3); n <- 200L
  df <- data.frame(a = factor(sample(c("A", "B"), n, TRUE)))
  list(df = df, tg = list(a = c(A = 0.5, B = 0.5)))
}

test_that("harvest rejects NaN / non-finite weight bounds (CR-D2)", {
  f <- .iv_fixture()
  expect_error(harvest(f$df, f$tg, min_weight = NaN), "min_weight")
  expect_error(harvest(f$df, f$tg, max_weight = NaN), "max_weight")
  expect_error(harvest(f$df, f$tg, min_weight = Inf), "min_weight")
  expect_error(harvest(f$df, f$tg, min_weight = NA_real_), "min_weight")
})

test_that("harvest still accepts +Inf max_weight (unbounded) (CR-D2)", {
  f <- .iv_fixture()
  expect_error(
    suppressWarnings(harvest(f$df, f$tg, max_weight = Inf, attach_weights = FALSE)),
    NA)   # no error
})

test_that("harvest rejects NA / non-positive iteration budgets (CR-D2)", {
  f <- .iv_fixture()
  expect_error(harvest(f$df, f$tg, max_iterations = NA_integer_), "max_iterations")
  expect_error(harvest(f$df, f$tg, max_iterations = -5L), "max_iterations")
  expect_error(harvest(f$df, f$tg, max_iterations = 0L), "max_iterations")
})
