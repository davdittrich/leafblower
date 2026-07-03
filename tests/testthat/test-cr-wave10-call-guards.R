# CR-D wave 10 direct-.Call hardening (subprocess-isolated so a regression crash
# surfaces as a test failure, not a whole-suite abort).
#   D13 (j7x8.13): capacity_penalty / alm_penalty passed as INTSXP (1L) must fall to
#     the -1.0 sentinel via a TYPEOF==REALSXP guard, NOT trigger a hard R accessor
#     crash. (Already guarded in current source; this locks the behavior in.)
#   D17 (j7x8.17): group_ids=list() (K=0) must produce a graceful R error, not reach
#     weights(n) outside the try and std::terminate.

# Build the 39-arg C_rk_calibrate argument list (mirrors R/harvest.R), overriding
# individual positions by index. Runs in a fresh subprocess.
.run_call <- function(overrides) {
  callr::r(function(ov) {
    library(leafblower)
    args <- list(
      list(as.integer(c(0L, 1L, 0L, 1L))),  # 1 group_ids (1 margin, 4 obs)
      2L,                                    # 2 cat_counts
      list(c(0.5, 0.5)),                     # 3 targets
      4L,                                    # 4 n_obs
      0.0, 5.0, "raking", 0L, 10L, NULL,     # 5..10
      -1.0, -1.0, 1e-6, 0L,                  # 11 capacity, 12 alm, 13 tol, 14 bounds_mode
      1L, 1.0, 1.0, 0.5,                     # 15..18 homotopy
      "round_robin", "fixed", 1.0, 1.0, 0.5, # 19..23
      0.001, 0.0, 0L, 1L, 0L,                # 24..28
      0L, 0L, 1.0, 0.3, 1.5, -1.0, 20L, 2L,  # 29..36
      0L, 1e-8, 0.0)                         # 37..39
    for (nm in names(ov)) args[[as.integer(nm)]] <- ov[[nm]]
    tryCatch({
      do.call(function(...) .Call("C_rk_calibrate", ..., PACKAGE = "leafblower"), args)
      "OK"
    }, error = function(e) paste0("ERR:", conditionMessage(e)))
  }, args = list(ov = overrides))
}

test_that("capacity_penalty / alm_penalty as INTSXP fall to sentinel, no crash (CR-D13)", {
  skip_if_not_installed("callr")
  # capacity_penalty = 1L (position 11), alm_penalty = 1L (position 12): TYPEOF guard
  # must reject non-REALSXP and use -1.0 sentinel -> raking runs to completion.
  expect_equal(.run_call(list("11" = 1L)), "OK")
  expect_equal(.run_call(list("12" = 1L)), "OK")
  # both at once
  expect_equal(.run_call(list("11" = 1L, "12" = 1L)), "OK")
})

test_that("group_ids=list() (K=0) errors gracefully, no crash (CR-D17)", {
  skip_if_not_installed("callr")
  msg <- .run_call(list("1" = list(), "2" = integer(0), "3" = list()))
  expect_match(msg, "at least one margin", fixed = TRUE)
})
