# tests/testthat/test-algo-selection.R
# Guard: prevent benchmark execution when sourced for testing
.BENCH_SOURCED <- TRUE

# Source the benchmark to load helper functions (no side effects with guard set)
# NOTE: benchmark file is built incrementally — only functions added so far are available.
# Use a simple heuristic: walk up from tests/testthat until we find DESCRIPTION
find_project_root <- function(start_path = getwd()) {
  path <- start_path
  while (path != dirname(path)) {  # Stop at filesystem root
    if (file.exists(file.path(path, "DESCRIPTION"))) {
      return(path)
    }
    path <- dirname(path)
  }
  stop("Project root (DESCRIPTION) not found")
}

bench_file <- file.path(find_project_root(), "benchmarks/algo_selection_benchmark.R")
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
  seeds2 <- sapply(c(-3.0, -4.0, -5.0, -6.0), function(x) bench_seed(5.0, x))
  expect_equal(length(unique(seeds2)), 4L)
})
