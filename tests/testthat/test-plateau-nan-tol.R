# 5ye4.17: a non-finite (NaN/Inf) plateau tol/pct must reject with the SAME
# canonical, bridge-neutral message as the Python bridge — not base-R's
# "missing value where TRUE/FALSE needed" (which fired before the range stop()).

test_that("plateau rule rejects a non-finite tol with the canonical message (5ye4.17)", {
  set.seed(1)
  df <- data.frame(a = factor(sample(c("x", "y"), 50, TRUE)))
  tg <- list(a = c(x = .5, y = .5))
  expect_error(
    harvest(df, tg, method = "raking", attach_weights = FALSE,
            convergence = list(tol = NaN, rule = "plateau")),
    "convergence tol must be a finite value in \\(0,1\\) for rule='plateau'")
  expect_error(
    harvest(df, tg, method = "raking", attach_weights = FALSE,
            convergence = list(tol = Inf, rule = "plateau")),
    "convergence tol must be a finite value in \\(0,1\\)")
})

test_that("plateau pct shorthand rejects a non-finite pct with the canonical message (5ye4.17)", {
  set.seed(1)
  df <- data.frame(a = factor(sample(c("x", "y"), 50, TRUE)))
  tg <- list(a = c(x = .5, y = .5))
  expect_error(
    harvest(df, tg, method = "raking", attach_weights = FALSE,
            convergence = list(pct = NaN)),
    "convergence pct must be a finite value in \\(0,1\\) for rule='plateau'")
})
