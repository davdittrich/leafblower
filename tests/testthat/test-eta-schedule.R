test_that("tang_dynamic eta reduces errRp vs fixed at matched budget", {
  skip_if_not_installed("arrow")
  fx <- test_path("fixtures/stepstone_small.parquet")
  tg <- test_path("fixtures/stepstone_small_targets.rds")
  skip_if(!file.exists(fx) || !file.exists(tg))
  data   <- arrow::read_parquet(fx)
  target <- readRDS(tg)

  common <- list(
    data = data, target = target,
    max_weight = 5, method = "ieppa",
    max_iterations = 500,
    convergence = list(absolute = 1e-12),
    homotopy_levels = 5,
    homotopy_start_factor = 10,
    homotopy_end_factor = 1,
    attach_weights = FALSE
  )
  probe_iters <- "50,100,200,500"
  tmp_fixed <- tempfile(fileext = ".csv")
  tmp_dyn   <- tempfile(fileext = ".csv")

  withr::with_envvar(
    c(LBW_TRAJECTORY_AT = probe_iters, LBW_TRAJECTORY_OUT = tmp_fixed),
    do.call(leafblower::harvest, c(common, list(eta_schedule = "fixed"))))

  withr::with_envvar(
    c(LBW_TRAJECTORY_AT = probe_iters, LBW_TRAJECTORY_OUT = tmp_dyn),
    do.call(leafblower::harvest, c(common, list(
      eta_schedule = "tang_dynamic",
      eta_start = 10, eta_end = 1,
      eta_schedule_power = 0.5))))

  fixed_err <- tail(utils::read.csv(tmp_fixed)$errRp, 1)
  dyn_err   <- tail(utils::read.csv(tmp_dyn)$errRp, 1)
  expect_lt(dyn_err, fixed_err)
})
