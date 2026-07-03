# CR-E13 (leafblower-5ye4.13): a long-format target with a duplicated (variable,
# level) row is ambiguous; R's parse_target (setNames) previously kept both under
# duplicate names. Both R and Python now raise on duplicate levels.

test_that("parse_target rejects duplicate (variable, level) rows (CR-E13)", {
  set.seed(1); n <- 100
  df <- data.frame(a = factor(sample(c("x", "y"), n, TRUE)))
  tdf <- data.frame(variable = c("a", "a", "a"),
                    level = c("x", "y", "x"),
                    proportion = c(0.5, 0.5, 0.3))
  expect_error(
    suppressWarnings(harvest(df, tdf, attach_weights = FALSE)),
    "duplicate level")
})

test_that("parse_target rejects NaN-level duplicates (CR-E13)", {
  set.seed(1); n <- 100
  df <- data.frame(a = factor(sample(c("x", "y"), n, TRUE)))
  tdf <- data.frame(variable = c("a", "a"), level = c(NA, NA),
                    proportion = c(0.5, 0.5))
  expect_error(
    suppressWarnings(harvest(df, tdf, attach_weights = FALSE)),
    "duplicate level")
})

test_that("parse_target accepts unique levels (CR-E13)", {
  set.seed(1); n <- 100
  df <- data.frame(a = factor(sample(c("x", "y"), n, TRUE)))
  tdf <- data.frame(variable = c("a", "a"), level = c("x", "y"),
                    proportion = c(0.5, 0.5))
  w <- suppressWarnings(harvest(df, tdf, attach_weights = FALSE))
  expect_length(w, n)
})
