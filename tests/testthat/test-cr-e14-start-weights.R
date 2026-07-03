# CR-E10b (leafblower-5ye4.14): R start_weights validation matches Python — reject
# negatives / non-finite; length-1 broadcasts to n; wrong length rejected.

test_that("normalize_start_weights rejects negatives and non-finite, broadcasts length-1 (CR-E10b)", {
  set.seed(1); n <- 100
  df <- data.frame(a = factor(sample(c("x", "y"), n, TRUE)))
  tg <- list(a = c(x = 0.5, y = 0.5))

  expect_error(
    suppressWarnings(harvest(df, tg, start_weights = c(rep(2, n - 1), -0.5), attach_weights = FALSE)),
    "negative")
  expect_error(
    suppressWarnings(harvest(df, tg, start_weights = c(rep(1, n - 1), NaN), attach_weights = FALSE)),
    "non-finite")
  expect_error(
    suppressWarnings(harvest(df, tg, start_weights = rep(1, 50), attach_weights = FALSE)),
    "length")
  expect_error(
    suppressWarnings(harvest(df, tg, start_weights = rep(1e-18, n), attach_weights = FALSE)),
    "positive value")
  expect_error(
    suppressWarnings(harvest(df, tg, start_weights = matrix(1, n, 1), attach_weights = FALSE)),
    "1-D")
  # length-1 broadcasts to n (parity with Python) — must NOT error, returns n weights
  w1 <- suppressWarnings(harvest(df, tg, start_weights = 3, attach_weights = FALSE))
  expect_length(w1, n)
  wn <- suppressWarnings(harvest(df, tg, start_weights = rep(1, n), attach_weights = FALSE))
  expect_length(wn, n)
})
