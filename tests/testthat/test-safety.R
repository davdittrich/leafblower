# tests/testthat/test-safety.R
library(leafblower)

.mini_harvest <- function(n = 20, method = "raking") {
  set.seed(1)
  df <- data.frame(
    sex  = sample(c("M", "F"), n, replace = TRUE),
    wt   = rep(1.0, n)
  )
  target <- list(sex = c(M = 0.5, F = 0.5))
  harvest(df, target = target, weight_column = "wt", method = method)
}

test_that("B4: solver dispatch returns without crashing for each method", {
  for (m in c("raking", "lbfgsb", "ieppa", "auto")) {
    expect_no_error(.mini_harvest(method = m))
  }
})

test_that("B7: C_rk_calibrate start_weights length mismatch produces error", {
  set.seed(2)
  n <- 30L
  df <- data.frame(
    sex = sample(c("M", "F"), n, replace = TRUE),
    wt  = rep(1.0, n)
  )
  targets <- list(sex = c(M = 0.5, F = 0.5))
  bad_sw  <- rep(1.0, n - 5L)   # wrong length: 25 instead of 30

  expect_error(
    .Call("C_rk_calibrate",
          df,                          # 1:  data
          targets,                     # 2:  target
          as.double(0),                # 3:  min_weight
          as.double(1e6),              # 4:  max_weight
          as.character("raking"),      # 5:  method
          as.integer(0L),              # 6:  verbose
          as.integer(50L),             # 7:  max_iterations
          bad_sw,                      # 8:  start_weights (wrong length — the probe)
          NULL,                        # 9:  capacity_penalty (NULL=auto; not under test here)
          as.double(1e-6),             # 10: tol_abs (legacy slot)
          as.integer(0L),              # 11: bounds_mode_int
          as.integer(1L),              # 12: homotopy_levels
          as.double(1.0),              # 13: homotopy_start_factor
          as.double(1.0),              # 14: homotopy_end_factor
          as.double(1.0),              # 15: homotopy_budget_p
          as.character("round_robin"), # 16: scheduler
          as.character("none"),        # 17: eta_schedule
          as.double(0.5),              # 18: eta_start
          as.double(0.1),              # 19: eta_end
          as.double(1.0),              # 20: eta_schedule_power
          as.double(1e-4),             # 21: conv pct_tol
          as.double(0.0),              # 22: conv absolute_tol
          as.integer(4L),              # 23: conv metric_int
          as.integer(0L),              # 24: conv rule_int
          as.integer(0L),              # 25: conv stop_when_int
          as.integer(0L),              # 26: sor enabled
          as.integer(0L),              # 27: sor auto
          as.double(1.0),              # 28: sor omega_init
          as.double(0.5),              # 29: sor omega_min
          as.double(1.0),              # 30: sor omega_fixed
          as.integer(0L),              # 31: sor burnin
          as.integer(0L)               # 32: accelerate_bool
    ),
    regexp = "start_weights length"
  )
})
