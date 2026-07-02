test_that("0ai8.1->eb79.23: logit tight-bounds infeasible reports INFEAS (was mislabeled STALL)", {
  library(leafblower)
  set.seed(42)
  # Perfectly correlated a=1 <-> b=x, a=2 <-> b=y (50 each). Target b=(x=0.8,y=0.2) with tight
  # bounds [0.999,1.001]: margin b=x needs 0.8*100=80 but only 50 obs are x, each capped near 1
  # => max mass ~50.05 < 80. STRUCTURALLY INFEASIBLE. This fixture previously exercised the
  # RK_ERR_STALL/stall_kl path (imprecise); eb79.23's pre-loop interval-sum check now correctly
  # reports RK_ERR_INFEAS => harvest() STOPS. The stall_kl convergence_reason mechanism itself
  # remains covered by test-stall-kind-greenkhorn.R.
  data <- data.frame(
    a = factor(c(rep("1", 50), rep("2", 50))),
    b = factor(c(rep("x", 50), rep("y", 50)))
  )
  target <- list(a = c("1" = 0.5, "2" = 0.5), b = c(x = 0.8, y = 0.2))
  expect_error(
    harvest(data, target, method = "logit", min_weight = 0.999, max_weight = 1.001,
            convergence = list(improvement = 1e-15), max_iterations = 5L, attach_weights = FALSE))
})
