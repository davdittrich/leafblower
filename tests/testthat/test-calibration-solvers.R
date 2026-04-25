test_that("T1a: new method names error cleanly (stubs)", {
  data <- data.frame(a = factor(c("1","2")))
  target <- list(a = c("1"=0.5, "2"=0.5))
  for (m in c("sinkhorn", "chebyshev", "greg", "grake")) {
    expect_error(
      leafblower::harvest(data, target, max_weight=3, method=m, attach_weights=FALSE),
      info = paste("method", m, "should error not crash")
    )
  }
})

test_that("T1b: convergence_used$objective and $minimized_metric present", {
  set.seed(1)
  data <- data.frame(a=factor(sample(c("1","2"),200,TRUE)))
  target <- list(a=c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=3, method="ieppa",
                           attach_weights=FALSE)
  r <- attr(w, "result")
  expect_true("objective" %in% names(r$convergence_used))
  expect_true("minimized_metric" %in% names(r$convergence_used))
  expect_true(is.finite(r$convergence_used$objective))
})
