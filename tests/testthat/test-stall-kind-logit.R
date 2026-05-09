test_that("0ai8.1: logit RK_ERR_STALL emits stall_kl convergence_reason", {
  library(leafblower)
  set.seed(42)
  # Perfectly correlated data: a=1 <-> b=x, a=2 <-> b=y.
  # Conflicting targets + tight bounds [0.999, 1.001] leave Newton with no
  # feasible direction: residual improvement < 0.1% → RK_ERR_STALL (stall_kind=2).
  data <- data.frame(
    a = factor(c(rep("1", 50), rep("2", 50))),
    b = factor(c(rep("x", 50), rep("y", 50)))
  )
  target <- list(a = c("1" = 0.5, "2" = 0.5), b = c(x = 0.8, y = 0.2))
  w <- harvest(
    data, target,
    method         = "logit",
    min_weight     = 0.999,
    max_weight     = 1.001,
    convergence    = list(improvement = 1e-15),
    max_iterations = 5L,
    attach_weights = FALSE
  )
  r <- attr(w, "result")
  expect_equal(r$status, 5L, label = "status must be RK_ERR_STALL=5")
  expect_equal(r$convergence_used$convergence_reason, "stall_kl")
})
