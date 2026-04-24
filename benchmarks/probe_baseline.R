data   <- arrow::read_parquet("tests/testthat/fixtures/stepstone_small.parquet")
target <- readRDS("tests/testthat/fixtures/stepstone_small_targets.rds")
out <- "benchmarks/baseline_trajectory_stepstone_small.csv"
Sys.setenv(LBW_TRAJECTORY_AT  = "1,10,50,100,200,500,1000,2000,3000",
           LBW_TRAJECTORY_OUT = out)
leafblower::harvest(
  data, target, max_weight = 5,
  method = "ieppa",
  max_iterations = 3000,
  convergence = list(absolute = 1e-12),
  attach_weights = FALSE
)
Sys.unsetenv(c("LBW_TRAJECTORY_AT", "LBW_TRAJECTORY_OUT"))
cat("Baseline trajectory written to", out, "\n")
