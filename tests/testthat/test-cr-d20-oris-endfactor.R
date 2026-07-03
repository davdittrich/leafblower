# CR-D20 (leafblower-j7x8.20): the oris preentry capacity check validates against the
# loosest homotopy final bound (max_weight * max(1, end_factor)) ONLY when homotopy is
# active (n_levels>1), plus a tightest-level min<=max guard -> no false-positive INFEAS
# when end_factor>1, without masking a genuine infeasibility.

test_that("oris homotopy end_factor>1 does not false-positive INFEAS (CR-D20)", {
  set.seed(5); n <- 300
  df <- data.frame(a  = factor(sample(c("x", "y"), n, TRUE)),
                   b  = factor(sample(c("p", "q"), n, TRUE)),
                   cc = factor(sample(c("m", "o"), n, TRUE)))
  tg <- list(a = c(x = .5, y = .5), b = c(p = .5, q = .5), cc = c(m = .5, o = .5))
  # max_weight=0.8 is infeasible at the BASE bound (sum_U=0.8n<n) but the homotopy
  # final level is 0.8*1.5=1.2 (sum_U=1.2n>=n) -> feasible. Must NOT INFEAS.
  r <- attr(suppressWarnings(harvest(df, tg, method = "oris", max_weight = 0.8,
              homotopy_levels = 3L, homotopy_start_factor = 1.0,
              homotopy_end_factor = 1.5)), "result")
  expect_equal(r$status, 0L)
  # Control 1: genuinely infeasible (final still 0.8n<n) still rejects.
  expect_error(
    suppressWarnings(harvest(df, tg, method = "oris", max_weight = 0.8,
      homotopy_levels = 3L, homotopy_start_factor = 1.0, homotopy_end_factor = 1.0,
      attach_weights = FALSE)))
  # Control 2: homotopy OFF (n_levels=1) never applies end_factor (driver forces
  # factor=1.0), so max_weight=0.8 is the effective bound and end_factor=1.5 must
  # NOT rescue it — the lift is gated on n_levels>1.
  expect_error(
    suppressWarnings(harvest(df, tg, method = "oris", max_weight = 0.8,
      homotopy_levels = 1L, homotopy_end_factor = 1.5, attach_weights = FALSE)))
})
