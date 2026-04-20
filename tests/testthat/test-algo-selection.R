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
