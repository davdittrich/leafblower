# tests/testthat/test-sraa-global.R
# T_sraa_global: SRAA-m global safeguard regression test (K=4 overlapping margins)

test_that("T_sraa_global: greenkhorn+SRAA max_err <= plain on K=4 overlapping-margin problem", {
  # K=4 designed to reproduce the basin-escape failure seen at K=9 (stepstone).
  # Old local safeguard: AA overshoots max_err-optimal basin on multi-margin problems.
  # New global safeguard: AA stays in or returns to max_err-optimal basin.
  # RED with old per-step local safeguard (err_AA <= err_plain).
  # GREEN after global best_err_seen safeguard + revert-to-best.
  set.seed(5); n <- 3000L
  df <- data.frame(
    a = factor(sample(letters[1:4], n, TRUE)),
    b = factor(sample(LETTERS[1:3], n, TRUE)),
    c = factor(sample(c("x", "y"),  n, TRUE)),
    d = factor(sample(c("M", "F"),  n, TRUE))
  )
  tgt <- list(
    a = setNames(c(0.4, 0.3, 0.2, 0.1), letters[1:4]),
    b = setNames(c(0.5, 0.3, 0.2),      LETTERS[1:3]),
    c = c(x = 0.6, y = 0.4),
    d = c(M = 0.45, F = 0.55)
  )
  r_aa    <- suppressWarnings(harvest(df, tgt, method = "greenkhorn", accelerate = TRUE,
                                       max_iterations = 500L, attach_weights = FALSE))
  r_plain <- suppressWarnings(harvest(df, tgt, method = "greenkhorn", accelerate = FALSE,
                                       max_iterations = 500L, attach_weights = FALSE))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label = sprintf(
      "SRAA K=4 overlapping-margin (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
})
