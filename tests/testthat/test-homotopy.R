test_that("P-A homotopy achieves >=30pct errRp reduction vs baseline (tight-clamp synthetic)", {
  set.seed(31415)
  n <- 5000
  K <- 5
  data <- data.frame(
    v1 = factor(sample(1:4, n, replace = TRUE)),
    v2 = factor(sample(1:3, n, replace = TRUE)),
    v3 = factor(sample(1:5, n, replace = TRUE)),
    v4 = factor(sample(1:2, n, replace = TRUE)),
    v5 = factor(sample(1:6, n, replace = TRUE))
  )
  target <- list(
    v1 = c("1"=0.1, "2"=0.4, "3"=0.4, "4"=0.1),
    v2 = c("1"=0.5, "2"=0.3, "3"=0.2),
    v3 = c("1"=0.1, "2"=0.1, "3"=0.4, "4"=0.3, "5"=0.1),
    v4 = c("1"=0.7, "2"=0.3),
    v5 = c("1"=0.05, "2"=0.05, "3"=0.5, "4"=0.2, "5"=0.15, "6"=0.05)
  )
  max_w <- 2  # tight clamp — many cells saturate, slow-rate regime

  tmp_base <- tempfile(fileext = ".csv")
  tmp_homo <- tempfile(fileext = ".csv")
  probe_iters <- "50,100,200,500"

  withr::with_envvar(
    c(LBW_TRAJECTORY_AT = probe_iters, LBW_TRAJECTORY_OUT = tmp_base),
    leafblower::harvest(data, target, max_weight = max_w,
                        method = "oris",
                        max_iterations = 500,
                        convergence = list(absolute = 1e-12),
                        attach_weights = FALSE))

  homo <- withr::with_envvar(
    c(LBW_TRAJECTORY_AT = probe_iters, LBW_TRAJECTORY_OUT = tmp_homo),
    leafblower::harvest(data, target, max_weight = max_w,
                        method = "oris",
                        max_iterations = 500,
                        convergence = list(absolute = 1e-12),
                        homotopy_levels = 5,
                        homotopy_start_factor = 10,
                        homotopy_end_factor = 1,
                        homotopy_budget_p = 0.5,
                        attach_weights = FALSE))

  base_err <- tail(utils::read.csv(tmp_base)$errRp, 1)
  homo_err <- tail(utils::read.csv(tmp_homo)$errRp, 1)

  # A2: homotopy must achieve >=30% errRp reduction at identical iteration budget
  expect_lt(homo_err, 0.7 * base_err)
})
