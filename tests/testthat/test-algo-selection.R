# tests/testthat/test-algo-selection.R
# Guard: prevent benchmark execution when sourced for testing
assign(".BENCH_SOURCED", TRUE, envir = .GlobalEnv)

# Source the benchmark to load helper functions (no side effects with guard set)
# NOTE: benchmark file is built incrementally — only functions added so far are available.
bench_file <- file.path(
  rprojroot::find_root(rprojroot::has_file("DESCRIPTION")),
  "benchmarks", "algo_selection_benchmark.R"
)
source(bench_file)

test_that("bench_seed is deterministic", {
  s1 <- bench_seed(5.0, -4.0)
  s2 <- bench_seed(5.0, -4.0)
  expect_identical(s1, s2)
})

test_that("bench_seed returns integer in 32-bit range", {
  s <- bench_seed(7.7, -3.0)
  expect_true(is.integer(s))
  expect_true(s > 0L && s < .Machine$integer.max)
})

test_that("bench_seed differs for different inputs", {
  seeds <- sapply(c(4.0, 5.0, 6.0, 7.0), function(x) bench_seed(x, -4.0))
  expect_equal(length(unique(seeds)), 4L)
  seeds2 <- sapply(c(-3.2, -4.5, -5.1, -5.9), function(x) bench_seed(5.0, x))
  expect_equal(length(unique(seeds2)), 4L)
})

test_that("make_bench_data returns correct structure", {
  set.seed(1L)
  out <- make_bench_data(n = 200L, K = 3L, cats_per_margin = 4L)
  expect_named(out, c("df", "targets"))
  expect_equal(nrow(out$df), 200L)
  expect_equal(ncol(out$df), 3L)
  expect_equal(length(out$targets), 3L)
})

test_that("make_bench_data targets sum to 1", {
  set.seed(2L)
  out <- make_bench_data(n = 100L, K = 2L, cats_per_margin = 5L)
  sums <- sapply(out$targets, sum)
  expect_true(all(abs(sums - 1.0) < 1e-12))
})

test_that("make_bench_data df columns are factors with correct levels", {
  set.seed(3L)
  out <- make_bench_data(n = 50L, K = 2L, cats_per_margin = 3L)
  expect_true(all(sapply(out$df, is.factor)))
  expect_true(all(sapply(out$df, nlevels) == 3L))
})

test_that("make_bench_data validates inputs", {
  expect_error(make_bench_data(n = 0L,  K = 2L, cats_per_margin = 3L))
  expect_error(make_bench_data(n = 10L, K = 0L, cats_per_margin = 3L))
  expect_error(make_bench_data(n = 10L, K = 2L, cats_per_margin = 1L))
})

test_that("make_bench_data df levels align with targets names", {
  set.seed(7L)
  out <- make_bench_data(n = 30L, K = 3L, cats_per_margin = 4L)
  for (k in seq_along(out$targets)) {
    expect_identical(names(out$targets[[k]]), levels(out$df[[k]]))
  }
})
