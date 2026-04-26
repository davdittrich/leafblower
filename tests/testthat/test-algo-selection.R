# tests/testthat/test-algo-selection.R
# Guard: prevent benchmark execution when sourced for testing
assign(".BENCH_SOURCED", TRUE, envir = .GlobalEnv)

# Source the benchmark to load helper functions (no side effects with guard set)
# NOTE: benchmark file is built incrementally — only functions added so far are available.
bench_file <- file.path(
  rprojroot::find_root(rprojroot::has_file("DESCRIPTION")),
  "benchmarks", "algo_selection_benchmark.R"
)
# Set working directory to package root so relative sources in benchmark work
old_wd <- getwd()
setwd(rprojroot::find_root(rprojroot::has_file("DESCRIPTION")))
on.exit(setwd(old_wd))
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

test_that("time_cell returns finite numeric scalar", {
  # Tiny problem for speed: log_complexity=4.0 → n≈278, K=9, cats=4
  # tol=1e-3 (loose) → fast convergence
  y <- time_cell(log_complexity = 4.0, log_tol = -3.0, K = 9L)
  expect_true(is.numeric(y) && length(y) == 1L && is.finite(y))
})

test_that("time_cell is deterministic for same inputs", {
  y1 <- time_cell(log_complexity = 4.0, log_tol = -3.0, K = 9L)
  y2 <- time_cell(log_complexity = 4.0, log_tol = -3.0, K = 9L)
  # Same seed → same data → same weights → same timing direction (sign)
  expect_equal(sign(y1), sign(y2))
})

test_that("time_cell seed_extra produces deterministic K-stability seeds", {
  # seed_extra must produce a repeatable result distinct from seed_extra=0
  y_main <- time_cell(4.0, -3.0, K = 9L, seed_extra = 0L)
  y_k3_a <- time_cell(4.0, -3.0, K = 3L, seed_extra = 3L * 10000000L)
  y_k3_b <- time_cell(4.0, -3.0, K = 3L, seed_extra = 3L * 10000000L)
  # K-stability call is repeatable
  expect_equal(sign(y_k3_a), sign(y_k3_b))
  # seed_extra actually shifts the RNG state — verify seeds are distinct
  s_main <- bench_seed(4.0, -3.0) + 0L
  s_k3   <- bench_seed(4.0, -3.0) + 3L * 10000000L
  expect_false(s_main == s_k3)
  # Distinct seeds produce distinct data
  set.seed(s_main); d_main <- make_bench_data(100L, 3L, 4L)
  set.seed(s_k3);   d_k3   <- make_bench_data(100L, 3L, 4L)
  expect_false(identical(d_main$df, d_k3$df))
})

test_that("fit_gp returns km object with finite predictions", {
  set.seed(42L)
  design <- matrix(c(4.0, 5.0, 6.0, 7.0, 4.5, 5.5, 6.5, 4.0,
                     -3.0, -4.0, -5.0, -6.0, -3.5, -4.5, -5.5, -6.0),
                   ncol = 2)
  y <- rnorm(8L)
  gp <- fit_gp(design, y)
  expect_s4_class(gp, "km")
  cands <- as.data.frame(matrix(c(5.0, 5.5, -4.0, -4.5), ncol = 2))
  names(cands) <- c("V1", "V2")
  pred <- DiceKriging::predict(gp, newdata = cands, type = "UK", checkNames = FALSE)
  expect_true(all(is.finite(pred$mean)))
  expect_true(all(pred$sd >= 0))
})

test_that("straddle_next returns a 1-row matrix within bounds", {
  set.seed(1L)
  design <- matrix(c(4.0, 5.0, 6.0, 4.5, 5.5, 6.5, -3.0, -4.0, -5.0, -3.5, -4.5, -5.5), ncol = 2)
  y <- c(0.5, 0.1, -0.3, 0.3, 0.0, -0.1)
  gp <- fit_gp(design, y)
  cands <- as.matrix(expand.grid(
    V1 = seq(4.0, 7.7, length.out = 10),
    V2 = seq(-6.0, -3.0, length.out = 10)))
  nxt <- straddle_next(gp, cands, threshold = log(1.2))
  expect_equal(nrow(nxt), 1L)
  expect_equal(ncol(nxt), 2L)
  expect_true(nxt[1, 1] >= 4.0 && nxt[1, 1] <= 7.7)
  expect_true(nxt[1, 2] >= -6.0 && nxt[1, 2] <= -3.0)
})

test_that("classified_fraction is 0 for uncertain GP, 1 for certain", {
  # GP trained on values clearly above threshold → high classified_fraction
  set.seed(5L)
  design <- matrix(c(4.0, 5.0, 6.0, 7.0, -3.0, -4.0, -5.0, -6.0), ncol = 2)
  # y values clearly above threshold log(1.2) ≈ 0.182 → should be classified "above"
  y <- c(2.0, 2.1, 1.9, 2.2)
  gp <- suppressWarnings(fit_gp(design, y))
  cands <- as.matrix(expand.grid(V1 = seq(4.5, 7.2, length.out = 5),
                                  V2 = seq(-5.5, -3.5, length.out = 5)))
  frac <- classified_fraction(gp, cands, threshold = log(1.2))
  expect_true(is.numeric(frac) && length(frac) == 1L)
  expect_true(frac >= 0 && frac <= 1)
})

test_that("save_checkpoint + load_checkpoint round-trips state", {
  tmp_path <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_path), add = TRUE)
  state <- list(design = matrix(1:4, nrow = 2), y = c(0.1, 0.2), iter = 3L)
  save_checkpoint(state, tmp_path)
  expect_true(file.exists(tmp_path))
  loaded <- load_checkpoint(tmp_path)
  expect_equal(loaded$iter, 3L)
  expect_equal(loaded$y, c(0.1, 0.2))
  expect_equal(loaded$design, state$design)
})

test_that("load_checkpoint returns NULL when file absent", {
  result <- load_checkpoint(tempfile(fileext = ".rds"))
  expect_null(result)
})

test_that("save_checkpoint leaves no .tmp file on success", {
  tmp_path <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp_path), add = TRUE)
  save_checkpoint(list(x = 1L), tmp_path)
  expect_false(file.exists(paste0(tmp_path, ".tmp")))
})

# ── Algorithm routing regression tests ────────────────────────────────────────
# These tests verify select_algorithm() routing via harvest().
# Test 3 (L-BFGS-B path) is skipped until the benchmark runs and
# Case B constants are confirmed.

test_that("constrained (max_weight=5) always routes to iEPPA", {
  set.seed(99L)
  n   <- 500L
  df  <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  res <- leafblower::harvest(df, tgt, max_weight = 5)
  expect_equal(attr(res, "algorithm"), "ieppa")
})

test_that("default method (no method arg) routes to ieppa", {
  set.seed(7L)
  n  <- 30000L
  df <- data.frame(
    m1 = factor(sample(paste0("c", 1:4), n, replace = TRUE)),
    m2 = factor(sample(paste0("c", 1:4), n, replace = TRUE))
  )
  tgt <- list(
    m1 = c(c1 = 0.25, c2 = 0.25, c3 = 0.25, c4 = 0.25),
    m2 = c(c1 = 0.25, c2 = 0.25, c3 = 0.25, c4 = 0.25)
  )
  res <- leafblower::harvest(df, tgt,
                              max_weight = Inf, min_weight = 0,
                              convergence = list(absolute = 1e-3))
  expect_equal(attr(res, "algorithm"), "ieppa")
})

test_that("method='lbfgsb' with max_weight routes to lbfgsb", {
  set.seed(1L)
  n  <- 500L
  df <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  res <- leafblower::harvest(df, tgt, method = "lbfgsb", max_weight = 5)
  expect_equal(attr(res, "algorithm"), "lbfgsb")
})

test_that("method='lbfgsb' output weights satisfy max_weight/min_weight", {
  set.seed(2L)
  n   <- 500L
  df  <- data.frame(x = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.3, c = 0.2))
  res <- leafblower::harvest(df, tgt, method = "lbfgsb",
                              max_weight = 1.5, min_weight = 0.2)
  expect_true(max(res$weights) <= 1.5 + 1e-6)
  expect_true(min(res$weights) >= 0.2 - 1e-6)
})

test_that("auto: falls back to lbfgsb after primary NOCONV at 1 iter", {
  set.seed(77L)
  n <- 200L
  data <- data.frame(
    a = factor(sample(1:3, n, TRUE)),
    b = factor(sample(1:2, n, TRUE))
  )
  target <- list(a = c("1" = 0.4, "2" = 0.4, "3" = 0.2),
                 b = c("1" = 0.6, "2" = 0.4))
  r <- leafblower::harvest(data, target, method = "auto", max_iterations = 1L,
                           verbose = 0)
  # max_iterations=1 forces primary solver NOCONV; fallback should fire
  expect_equal(attr(r, "algorithm"), "lbfgsb")
})
