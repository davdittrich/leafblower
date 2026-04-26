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

test_that("T-routing: AUTO uses raking when incompressible (M_cell/n > 0.9)", {
  # 200 obs each in its own group → M_cell = 200, ratio = 1.0 (fully incompressible).
  # AUTO must select raking (not iEPPA) because there is no cell compression benefit.
  set.seed(99)
  n2 <- 200
  # Each obs is its own group for margin 'a' → M_cell = n2
  data <- data.frame(
    a = factor(seq_len(n2)),
    b = factor(sample(c("1", "2"), n2, replace = TRUE))
  )
  target <- list(
    a = setNames(rep(1 / n2, n2), as.character(seq_len(n2))),
    b = c("1" = 0.5, "2" = 0.5)
  )
  w <- leafblower::harvest(data, target, max_weight = 5, method = "auto",
                           max_iterations = 200, attach_weights = FALSE)
  r <- attr(w, "result")
  expect_equal(r$status, 0L,
               label = "AUTO must converge even with incompressible data")
  expect_lt(r$max_error, 0.1)
  # Routing check: RK_ALG_RAKING == 3
  expect_equal(r$algorithm_used, 3L,
               label = "AUTO must select raking when M_cell/n > 0.9")
})

test_that("A5: raking cell-table: correct + speedup vs obs-level reference", {
  skip_on_cran()
  skip_if_not_installed("arrow")
  skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  ref_path <- test_path("fixtures/raking_obs_reference_stepstone.rds")
  skip_if(!file.exists(ref_path), "obs-level raking reference fixture not generated")
  ref <- readRDS(ref_path)

  library(arrow); library(jsonlite)
  data <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  tgt_raw <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
  target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })

  t0 <- proc.time()["elapsed"]
  suppressWarnings(
    w_cell <- leafblower::harvest(data, target, method="raking",
                                  max_weight=5, max_iterations=500, attach_weights=FALSE)
  )
  elapsed_cell <- proc.time()["elapsed"] - t0

  r_cell <- attr(w_cell, "result")

  # Convergence
  expect_equal(r_cell$status, 0L, label="raking cell-table must converge")

  # Correctness: max_err within 1e-8 of obs-level reference (homogeneous d_i)
  expect_lt(abs(r_cell$max_error - ref$max_error), 1e-8,
            label="cell-table max_err must match obs-level to 1e-8")

  # Speedup: < 5% of obs-level elapsed (spec A5: elapsed_cell < elapsed_obs * 0.05)
  expect_lt(elapsed_cell, ref$elapsed_obs * 0.05,
            label=sprintf("cell-table must be < 5%% of obs-level elapsed (obs=%.1fs)", ref$elapsed_obs))

  # Hard bounds: all weights in [min_weight, max_weight]
  expect_true(all(w_cell >= 1/5 - 1e-10 & w_cell <= 5 + 1e-10),
              label="all weights within [1/5, 5] bounds")
})

test_that("R-bounds: raking respects min_weight/max_weight exactly after fix", {
  # Tight bounds + skewed targets force clamping; normalization after clamp would violate bounds.
  set.seed(7)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2","3","4","5"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
  target <- list(
    a = c("1"=0.5, "2"=0.2, "3"=0.15, "4"=0.1, "5"=0.05),
    b = c("1"=0.8, "2"=0.2)
  )
  out <- leafblower::harvest(data, target, min_weight=0.5, max_weight=1.5,
                             method="raking", max_iterations=500, attach_weights=TRUE)
  w <- out$weights
  expect_true(all(w >= 0.5 - 1e-10),
              info=sprintf("min weight %.6f < 0.5", min(w)))
  expect_true(all(w <= 1.5 + 1e-10),
              info=sprintf("max weight %.6f > 1.5", max(w)))
})

test_that("T-overflow: AUTO routing + algorithm_used populated", {
  set.seed(1)
  n <- 1e5L
  data <- data.frame(
    a = factor(sample(1:100, n, replace=TRUE)),
    b = factor(sample(1:2, n, replace=TRUE))
  )
  target <- list(
    a = setNames(rep(0.01, 100), as.character(1:100)),
    b = c("1"=0.5, "2"=0.5)
  )
  w <- leafblower::harvest(data, target, max_weight=5, method="auto", attach_weights=FALSE)
  r <- attr(w, "result")
  expect_equal(r$status, 0L)
  expect_lt(r$max_error, 0.01)
  # algorithm_used: AUTO with M_cell << n should route to ieppa (algorithm=1)
  expect_equal(r$algorithm_used, 1L,
               info="AUTO must select ieppa when M_cell/n << 0.9")
})
