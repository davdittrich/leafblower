test_that("LBW_TRAJECTORY_AT env var produces probe file", {
  skip_if_not_installed("arrow")
  fx <- test_path("fixtures/stepstone_small.parquet")
  tg <- test_path("fixtures/stepstone_small_targets.rds")
  skip_if(!file.exists(fx) || !file.exists(tg))
  data   <- arrow::read_parquet(fx)
  target <- readRDS(tg)
  out <- tempfile(fileext = ".csv")
  withr::with_envvar(
    c(LBW_TRAJECTORY_AT  = "1,10,50,100,200,500",
      LBW_TRAJECTORY_OUT = out),
    {
      leafblower::harvest(
        data, target, max_weight = 5,
        method = "ieppa",
        max_iterations = 500,
        convergence = list(absolute = 1e-12),
        attach_weights = FALSE
      )
    }
  )
  expect_true(file.exists(out))
  probes <- utils::read.csv(out)
  expect_named(probes, c("iter", "errRp"))
  expect_equal(probes$iter, c(1L, 10L, 50L, 100L, 200L, 500L))
  expect_length(probes$errRp, 6)
  expect_true(all(is.finite(probes$errRp) & probes$errRp > 0))
})
