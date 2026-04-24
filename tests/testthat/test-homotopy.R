test_that("P-A homotopy reduces errRp on stepstone-small", {
  skip_if_not_installed("arrow")
  fx <- test_path("fixtures/stepstone_small.parquet")
  tg <- test_path("fixtures/stepstone_small_targets.rds")
  skip_if(!file.exists(fx) || !file.exists(tg))
  data   <- arrow::read_parquet(fx)
  target <- readRDS(tg)

  probe_iters <- "1,10,50,100,200,500"
  tmp_base <- tempfile(fileext = ".csv")
  tmp_homo <- tempfile(fileext = ".csv")

  withr::with_envvar(
    c(LBW_TRAJECTORY_AT = probe_iters, LBW_TRAJECTORY_OUT = tmp_base),
    leafblower::harvest(data, target, max_weight = 5,
                        method = "ieppa",
                        max_iterations = 500,
                        convergence = list(absolute = 1e-12),
                        attach_weights = FALSE))

  homo <- withr::with_envvar(
    c(LBW_TRAJECTORY_AT = probe_iters, LBW_TRAJECTORY_OUT = tmp_homo),
    leafblower::harvest(data, target, max_weight = 5,
                        method = "ieppa",
                        max_iterations = 500,
                        convergence = list(absolute = 1e-12),
                        homotopy_levels = 5,
                        homotopy_start_factor = 10,
                        homotopy_end_factor = 1,
                        homotopy_budget_p = 0.5,
                        attach_weights = FALSE))

  # weights must respect final target max_weight
  expect_true(max(as.numeric(homo)) <= 5 + 1e-10)

  # homotopy should reach lower terminal errRp (A2: >=30% reduction)
  base_err <- tail(utils::read.csv(tmp_base)$errRp, 1)
  homo_err <- tail(utils::read.csv(tmp_homo)$errRp, 1)
  expect_lt(homo_err, 0.7 * base_err)
})
