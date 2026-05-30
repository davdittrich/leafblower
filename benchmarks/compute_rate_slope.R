#!/usr/bin/env Rscript
# benchmarks/compute_rate_slope.R
#
# WU-6: Rate-slope probe on stepstone-small (ABE trajectory).
# Gate: slope <= -0.75 (beats O(t^-0.5) baseline).

suppressPackageStartupMessages(library(leafblower))

data_s   <- arrow::read_parquet("tests/testthat/fixtures/stepstone_small.parquet")
target_s <- readRDS("tests/testthat/fixtures/stepstone_small_targets.rds")
out_abe  <- "benchmarks/ABE_trajectory_stepstone_small.csv"

cat("=== Rate-slope probe (stepstone-small, ABE config) ===\n")
cat(sprintf("n = %s rows, %d margins\n",
    format(nrow(data_s), big.mark = ","), length(target_s)))

Sys.setenv(
  LBW_TRAJECTORY_AT  = paste(c(1, 10, 50, 100, 200, 500, 1000, 2000, 3000), collapse = ","),
  LBW_TRAJECTORY_OUT = out_abe
)

suppressWarnings(
  leafblower::harvest(
    data_s, target_s,
    max_weight          = 5,
    method              = "oris",
    max_iterations      = 3000,
    convergence         = list(absolute = 1e-12),
    homotopy_levels     = 5L,
    homotopy_start_factor = 10.0,
    homotopy_end_factor = 1.0,
    scheduler           = "greedy",
    eta_schedule        = "tang_dynamic",
    eta_start           = 20.0,
    eta_end             = 1.0,
    eta_schedule_power  = 0.5,
    attach_weights      = FALSE
  )
)
Sys.unsetenv(c("LBW_TRAJECTORY_AT", "LBW_TRAJECTORY_OUT"))

if (!file.exists(out_abe)) {
  stop("Trajectory file not written — LBW_TRAJECTORY_AT env var may not be supported yet")
}

probe <- utils::read.csv(out_abe)
cat("\nTrajectory:\n")
print(probe)

fit   <- stats::lm(log(errRp) ~ log(iter), data = probe)
slope <- stats::coef(fit)["log(iter)"]
cat(sprintf("\nRate slope: %.3f (gate: <= -0.75)\n", slope))

saveRDS(list(slope = slope, fit = fit, probe = probe),
        "benchmarks/final_rate_slope.rds")

if (slope <= -0.75) {
  cat("GATE: PASS\n")
} else {
  cat(sprintf("GATE: FAIL (slope=%.3f > -0.75)\n", slope))
}

stopifnot(slope <= -0.75)
cat("Done.\n")
