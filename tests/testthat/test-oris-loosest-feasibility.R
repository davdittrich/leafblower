# j7x8.23: ORIS preentry total-capacity check validates against the LOOSEST homotopy
# level (mw_save * max(start_factor, end_factor)). The old max(1.0, end_factor) form
# floored the bound to base and ignored start_factor, so a start>end schedule with
# max_weight<1 (whose loose early levels reach Σw=n via homotopy_break) was
# false-rejected as INFEAS. A problem infeasible even at the loosest level still is.

test_that("start>end homotopy is NOT false-rejected when the loosest level is feasible (j7x8.23)", {
  set.seed(5); n <- 2000L
  df <- data.frame(a = factor(sample(c("p", "q", "r"), n, TRUE, prob = c(.5, .4, .1))))
  tg <- list(a = c(p = 1/3, q = 1/3, r = 1/3))
  # max_weight=0.6, start=2.0, end=1.0: loosest cap = 0.6*2.0 = 1.2 >= 1 (total mass n
  # reachable at the loose early levels), tightest = 0.6. The old floor-to-1 bound
  # (0.6*max(1,1)=0.6 < 1) rejected this as INFEAS; the loosest-level check does not.
  w <- suppressWarnings(harvest(df, tg, method = "oris", max_weight = 0.6, min_weight = 0.1,
                                homotopy_levels = 4L, homotopy_start_factor = 2.0,
                                homotopy_end_factor = 1.0, max_iterations = 300L,
                                attach_weights = FALSE))
  r <- attr(w, "result")
  expect_false(r$status == 2L)   # not pre-rejected INFEAS (loosest bound reaches Σw=n)
})

test_that("homotopy IS rejected when infeasible even at the loosest level (j7x8.23)", {
  set.seed(5); n <- 2000L
  df <- data.frame(a = factor(sample(c("p", "q", "r"), n, TRUE, prob = c(.5, .4, .1))))
  tg <- list(a = c(p = 1/3, q = 1/3, r = 1/3))
  # max_weight=0.4, start=2.0, end=1.0: loosest cap = 0.4*2.0 = 0.8 < 1 → Σw=n
  # unreachable at EVERY scheduled level → genuine INFEAS (harvest stops).
  expect_error(
    harvest(df, tg, method = "oris", max_weight = 0.4, min_weight = 0.1,
            homotopy_levels = 4L, homotopy_start_factor = 2.0, homotopy_end_factor = 1.0,
            max_iterations = 300L, attach_weights = FALSE))
})
