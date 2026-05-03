test_that("A5: best_error <= max_error for iEPPA", {
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1", "2"), n, replace = TRUE)),
    b = factor(sample(c("1", "2"), n, replace = TRUE))
  )
  target <- list(a = c("1" = 0.5, "2" = 0.5), b = c("1" = 0.5, "2" = 0.5))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 200,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_lte(result$best_error, result$max_error)
  expect_true(is.finite(result$best_error))
  expect_true(result$best_iter > 0)
})

test_that("A5: best_error <= max_error for raking", {
  set.seed(43)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1", "2"), n, replace = TRUE)),
    b = factor(sample(c("1", "2"), n, replace = TRUE))
  )
  target <- list(a = c("1" = 0.5, "2" = 0.5), b = c("1" = 0.5, "2" = 0.5))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "raking",
                           max_iterations = 200,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_lte(result$best_error, result$max_error)
})

test_that("A5: best_weights is obs-level, length n, sum normalized to n", {
  set.seed(44)
  n <- 1000
  data <- data.frame(a = factor(sample(c("1", "2"), n, replace = TRUE)))
  target <- list(a = c("1" = 0.5, "2" = 0.5))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_length(result$best_weights, n)
  expect_equal(sum(result$best_weights), n, tolerance = 1e-6)
})

test_that("A6: stepstone best-iterate within 5% of reference", {
  skip_on_cran()
  ref_path <- test_path("fixtures/stepstone_best_error_ref.rds")
  skip_if(!file.exists(ref_path))
  fx <- "benchmarks/stepstone_fulldata_bench_data.parquet"
  tg <- "benchmarks/stepstone_fulldata_bench_targets.json"
  skip_if(!file.exists(fx) || !file.exists(tg))
  data <- arrow::read_parquet(fx)
  tgt_raw <- jsonlite::fromJSON(tg)
  target <- lapply(tgt_raw, function(x) {
    v <- unlist(x)
    v / sum(v)
  })
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 3000,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  ref <- readRDS(ref_path)
  expect_lte(result$best_error, ref * 1.05)
})

test_that("z8wx: best_weights sum=n and best_error<=max_error with homotopy_levels=3", {
  set.seed(77)
  n <- 800
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3","4"), n, replace = TRUE))
  )
  target <- list(
    a = c("1"=0.4, "2"=0.4, "3"=0.2),
    b = c("1"=0.25, "2"=0.25, "3"=0.25, "4"=0.25)
  )
  w <- leafblower::harvest(data, target, max_weight = 2, method = "ieppa",
                           max_iterations = 300,
                           homotopy_levels = 3L,
                           homotopy_start_factor = 4.0,
                           homotopy_end_factor = 1.0,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(sum(result$best_weights), n, tolerance = 1e-6)
  expect_lte(result$best_error, result$max_error)
  expect_true(all(result$best_weights >= 0))
})

test_that("qbsf: best_iter is positive (cumulative counter confirmed)", {
  set.seed(88)
  n <- 500
  data <- data.frame(a = factor(sample(c("1","2","3"), n, replace = TRUE)))
  target <- list(a = c("1"=0.4, "2"=0.4, "3"=0.2))
  w <- leafblower::harvest(data, target, max_weight = 2, method = "ieppa",
                           max_iterations = 100,
                           convergence = list(absolute = 1e-10),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_gte(result$best_iter, 1L)
  expect_true(is.numeric(result$best_iter) || is.integer(result$best_iter))
})

test_that("g4oj: best_iter tracks active metric (kl vs max_err differ)", {
  # Use raking with max_iterations=3 (no convergence) so both solvers stop mid-run,
  # ensuring KL != errRp at the best snapshot.  iEPPA converges in <5 iters on this
  # dataset, collapsing both metrics to float-zero and making the assertion vacuous.
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.4,"2"=0.4,"3"=0.2), b = c("1"=0.6,"2"=0.4))
  suppressWarnings({
    w_kl <- leafblower::harvest(data, target, max_weight = 5, method = "raking",
                                max_iterations = 3,
                                convergence = list(metric="kl", rule="improvement", tol=1e-15),
                                attach_weights = FALSE)
    w_me <- leafblower::harvest(data, target, max_weight = 5, method = "raking",
                                max_iterations = 3,
                                convergence = list(metric="max_err", rule="improvement", tol=1e-15),
                                attach_weights = FALSE)
  })
  r_kl <- attr(w_kl, "result")
  r_me <- attr(w_me, "result")
  # KL run: best_error must equal the KL value (not errRp)
  expect_true(is.finite(r_kl$best_error))
  expect_gt(r_kl$best_error, 0)
  expect_equal(r_kl$best_error, r_kl$kl,
               info = "KL run: best_error should equal kl metric at best snapshot")
  # max_err run: best_error must equal max_error (errRp at best snapshot)
  expect_lte(r_me$best_error, r_me$max_error)
  # The two runs must produce strictly unequal best_error values (KL << errRp in absolute terms)
  expect_true(r_kl$best_error != r_me$best_error,
              info = "KL best_error != max_err best_error (different metrics tracked)")
})
