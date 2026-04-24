test_that("A5: best_error <= max_error for iEPPA", {
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1", "2"), n, replace = TRUE)),
    b = factor(sample(c("1", "2"), n, replace = TRUE))
  )
  target <- list(a = c("1" = 0.5, "2" = 0.5), b = c("1" = 0.5, "2" = 0.5))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 200, attach_weights = FALSE)
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
                           max_iterations = 200, attach_weights = FALSE)
  result <- attr(w, "result")
  expect_lte(result$best_error, result$max_error)
})

test_that("A5: best_weights is obs-level, length n, sum normalized to n", {
  set.seed(44)
  n <- 1000
  data <- data.frame(a = factor(sample(c("1", "2"), n, replace = TRUE)))
  target <- list(a = c("1" = 0.5, "2" = 0.5))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
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
                           max_iterations = 3000, attach_weights = FALSE)
  result <- attr(w, "result")
  ref <- readRDS(ref_path)
  expect_lte(result$best_error, ref * 1.05)
})
