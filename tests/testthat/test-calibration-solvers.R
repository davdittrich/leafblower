test_that("T1a: new method names error cleanly (stubs)", {
  data <- data.frame(a = factor(c("1","2")))
  target <- list(a = c("1"=0.5, "2"=0.5))
  for (m in c("sinkhorn", "chebyshev", "greg", "grake")) {
    expect_error(
      leafblower::harvest(data, target, max_weight=3, method=m, attach_weights=FALSE),
      regexp = "not yet implemented",
      info = paste("method", m, "should error with 'not yet implemented'")
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

test_that("T3: ieppa default convergence is kl+improvement", {
  set.seed(42)
  data <- data.frame(
    a = factor(sample(c("1","2","3"), 500, replace=TRUE)),
    b = factor(sample(c("1","2"), 500, replace=TRUE))
  )
  target <- list(a=c("1"=1/3,"2"=1/3,"3"=1/3), b=c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=10, method="ieppa",
                           max_iterations=500, attach_weights=FALSE)
  r <- attr(w, "result")
  expect_equal(r$convergence_used$metric, "kl",
               info="ieppa default metric must be 'kl' after spec change")
  expect_equal(r$convergence_used$rule, "improvement")
})

test_that("A1: sinkhorn KL <= ieppa KL at best_iter", {
  skip_on_cran()
  skip_if_not_installed("arrow")
  skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  ref_path <- test_path("fixtures/ieppa_kl_reference_stepstone.rds")
  skip_if(!file.exists(ref_path), "ieppa KL reference fixture not generated yet")
  ref <- readRDS(ref_path)

  # No unconditional skip — test FAILS with "sinkhorn not yet implemented"
  # until Plan C implements method="sinkhorn". That failure IS the RED state.
  data <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  tgt_raw <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
  target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })

  # Plan A: this expect_no_error call FAILS (sinkhorn is a stub returning RK_ERR_BADARG)
  # Plan C GREEN: sinkhorn is implemented so expect_no_error passes
  expect_no_error(
    w_s <- leafblower::harvest(data, target, method="sinkhorn",
                                max_weight=5, max_iterations=3000,
                                attach_weights=FALSE),
    message="Plan C must implement method='sinkhorn' to reach this point"
  )
  r_s <- attr(w_s, "result")
  expect_lte(r_s$convergence_used$objective, ref$kl_at_best_iter,
             label="sinkhorn KL <= ieppa best_iter KL")
  expect_equal(r_s$status, 0L, label="sinkhorn must converge")
  expect_equal(r_s$convergence_used$fired_at_iter, r_s$iterations,
               label="A7: sinkhorn best_iter == last_iter (monotone)")
})
