context("ieppa Jacobi log-path sweep")

test_that("Jacobi log-path matches GS within tolerance (K=3)", {
  skip_if_not_installed("withr")
  set.seed(42L)
  n <- 2000L
  df <- data.frame(
    a = factor(sample(letters[1:5], n, TRUE)),
    b = factor(sample(letters[1:4], n, TRUE)),
    c = factor(sample(letters[1:3], n, TRUE))
  )
  tgt <- list(
    a = setNames(rep(1/5, 5L), letters[1:5]),
    b = setNames(rep(1/4, 4L), letters[1:4]),
    c = setNames(rep(1/3, 3L), letters[1:3])
  )
  withr::with_envvar(c(LBW_IEPPA_FORCE_PATH = "log"), {
    r_gs <- suppressWarnings(harvest(df, tgt, method = "ieppa",
                                     convergence = list(absolute = 1e-4),
                                     jacobi_sweep = FALSE))
    r_ja <- suppressWarnings(harvest(df, tgt, method = "ieppa",
                                     convergence = list(absolute = 1e-4),
                                     jacobi_sweep = TRUE))
  })
  rg <- attr(r_gs, "result"); rj <- attr(r_ja, "result")
  expect_true(rg$status %in% c(0L, 4L, 5L))
  expect_true(rj$status %in% c(0L, 4L, 5L))
  expect_lt(abs(rg$max_error - rj$max_error), 2e-4,
    label = sprintf("GS max_err=%.4e vs Jacobi max_err=%.4e", rg$max_error, rj$max_error))
  expect_equal(sum(r_ja$weights), n, tolerance = 1e-6)
})

test_that("Jacobi log-path matches GS within tolerance (K=4)", {
  skip_if_not_installed("withr")
  set.seed(43L)
  n <- 2500L
  df <- data.frame(
    a = factor(sample(letters[1:5], n, TRUE)),
    b = factor(sample(letters[1:4], n, TRUE)),
    c = factor(sample(letters[1:3], n, TRUE)),
    d = factor(sample(letters[1:6], n, TRUE))
  )
  tgt <- list(
    a = setNames(rep(1/5, 5L), letters[1:5]),
    b = setNames(rep(1/4, 4L), letters[1:4]),
    c = setNames(rep(1/3, 3L), letters[1:3]),
    d = setNames(rep(1/6, 6L), letters[1:6])
  )
  withr::with_envvar(c(LBW_IEPPA_FORCE_PATH = "log"), {
    r_gs <- suppressWarnings(harvest(df, tgt, method = "ieppa",
                                     convergence = list(absolute = 1e-4),
                                     jacobi_sweep = FALSE))
    r_ja <- suppressWarnings(harvest(df, tgt, method = "ieppa",
                                     convergence = list(absolute = 1e-4),
                                     jacobi_sweep = TRUE))
  })
  rg <- attr(r_gs, "result"); rj <- attr(r_ja, "result")
  expect_true(rg$status %in% c(0L, 4L, 5L))
  expect_true(rj$status %in% c(0L, 4L, 5L))
  expect_lt(abs(rg$max_error - rj$max_error), 2e-4,
    label = sprintf("GS max_err=%.4e vs Jacobi max_err=%.4e", rg$max_error, rj$max_error))
  expect_equal(sum(r_ja$weights), n, tolerance = 1e-6)
})
