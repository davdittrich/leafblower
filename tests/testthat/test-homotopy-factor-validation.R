# j7x8.21: reject non-finite / non-positive homotopy factors at the .Call boundary.
# A NaN/Inf factor flows into the geometric schedule k_start*pow(k_end/k_start,frac)
# → NaN current_max_weight, and a NaN start_factor silently disables the CR-D20
# tightest-level feasibility guard.

test_that("harvest rejects non-finite / non-positive homotopy factors (j7x8.21)", {
  set.seed(1)
  df <- data.frame(a = factor(sample(c("x", "y"), 100, TRUE)))
  tg <- list(a = c(x = .5, y = .5))
  run <- function(...) harvest(df, tg, method = "oris", homotopy_levels = 2L,
                               attach_weights = FALSE, ...)
  for (bad in list(NaN, Inf, -Inf, 0, -1)) {
    expect_error(run(homotopy_start_factor = bad), "finite and > 0")
    expect_error(run(homotopy_end_factor   = bad), "finite and > 0")
  }
  # a valid tightening schedule still runs (guard not over-broad)
  expect_silent(suppressWarnings(
    run(homotopy_start_factor = 2.0, homotopy_end_factor = 1.0)))
})
