# CR-D7 (leafblower-j7x8.7): a direct .Call to C_rk_calibrate with n_obs <= 0
# (e.g. -1 or NA_integer_) must produce a graceful R error, NOT a std::terminate
# / R session abort. The std::vector<double> weights(n) allocation sits outside the
# solver try-block, so a negative/NA n would otherwise throw length_error/bad_alloc
# uncaught. Unreachable via harvest() (it enforces nrow >= 1) — direct-.Call only.
# Run in a subprocess so a regression (crash) surfaces as a distinct test failure
# rather than aborting the whole test run.

test_that("direct .Call with n_obs=-1 errors gracefully, no crash (CR-D7)", {
  skip_if_not_installed("callr")

  run_bad_nobs <- function(bad_n) {
    callr::r(function(nbad) {
      library(leafblower)
      # 39-arg C_rk_calibrate call mirroring R/harvest.R; only n_obs is invalid.
      args <- list(
        list(as.integer(c(0L, 1L, 0L, 1L))),  # 1 group_ids: one margin, 4 obs
        2L,                                    # 2 cat_counts
        list(c(0.5, 0.5)),                     # 3 targets
        nbad,                                  # 4 n_obs  <-- invalid
        0.0, 5.0, "raking", 0L, 10L, NULL,     # 5..10 min/max/method/verbose/maxit/sw
        -1.0, -1.0, 1e-6, 0L,                  # 11..14 capacity/alm/tol/bounds_mode
        1L, 1.0, 1.0, 0.5,                     # 15..18 homotopy levels/start/end/budget
        "round_robin", "fixed", 1.0, 1.0, 0.5, # 19..23 scheduler/eta*
        0.001, 0.0, 0L, 1L, 0L,                # 24..28 pct/abs/metric/rule/stop_when
        0L, 0L, 1.0, 0.3, 1.5, -1.0, 20L, 2L,  # 29..36 sor enabled/auto/omegas/burnin/mode
        0L, 1e-8, 0.0)                         # 37..39 accelerate/tsvd/ridge
      tryCatch({
        do.call(function(...) .Call("C_rk_calibrate", ..., PACKAGE = "leafblower"), args)
        "NO_ERROR"
      }, error = function(e) paste0("ERR:", conditionMessage(e)))
    }, args = list(nbad = bad_n))
  }

  msg_neg <- run_bad_nobs(-1L)
  expect_match(msg_neg, "n_obs must be a positive integer", fixed = TRUE)

  msg_na <- run_bad_nobs(NA_integer_)
  expect_match(msg_na, "n_obs must be a positive integer", fixed = TRUE)
})
