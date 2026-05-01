context("ieppa SRAA-m log-path extension")

## Test 1 — SRAA activates on forced log path ----------------------------

test_that("SRAA-m activates on log path (LBW_IEPPA_FORCE_PATH=log)", {
  skip_if_not_installed("withr")
  set.seed(2026L)
  n <- 20000L
  df  <- data.frame(
    a = factor(sample(4L, n, replace = TRUE)),
    b = factor(sample(5L, n, replace = TRUE))
  )
  tgt <- list(
    a = setNames(rep(0.25, 4L), as.character(1:4)),
    b = setNames(rep(0.2,  5L), as.character(1:5))
  )
  withr::with_envvar(c(LBW_IEPPA_FORCE_PATH = "log"), {
    r <- suppressWarnings(
      harvest(df, tgt, method = "ieppa",
              max_iterations = 200L, accelerate = TRUE)
    )
  })
  res <- attr(r, "result")
  expect_true(is.finite(res$max_error),
              label = "log-path SRAA: finite max_error")
  expect_equal(sum(r$weights), n, tolerance = 1e-6,
               label = "log-path SRAA: sum(weights) == n")
  expect_gte(res$aa_accepted_count, 1L,
             label = sprintf("SRAA must accept >= 1 log-path step (got %d)",
                             res$aa_accepted_count))
})

## Test 2 — Log-path SRAA validity (stepstone_small) ---------------------

test_that("SRAA-m log path produces valid weights on stepstone_small", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("withr")
  pq  <- testthat::test_path("fixtures/stepstone_small.parquet")
  rds <- testthat::test_path("fixtures/stepstone_small_targets.rds")
  skip_if_not(file.exists(pq))
  skip_if_not(file.exists(rds))
  ss  <- arrow::read_parquet(pq)
  tgt <- readRDS(rds)
  for (nm in names(tgt)) ss[[nm]] <- factor(ss[[nm]])
  withr::with_envvar(c(LBW_IEPPA_FORCE_PATH = "log"), {
    r <- suppressWarnings(
      harvest(ss, tgt, method = "ieppa",
              max_iterations = 100L, accelerate = TRUE)
    )
  })
  res <- attr(r, "result")
  expect_true(is.finite(res$max_error),
              label = "stepstone log SRAA: finite max_error")
  expect_equal(sum(r$weights), nrow(ss), tolerance = 1e-6,
               label = "stepstone log SRAA: sum == n")
  expect_true(res$status %in% c(0L, 1L, 4L, 5L),
              label = sprintf("valid status %d", res$status))
})

## Test 3 — Linear-to-log fallback with accelerate=TRUE ------------------

test_that("ieppa SRAA: linear-to-log fallback with accelerate=TRUE stays valid", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("withr")
  pq  <- testthat::test_path("fixtures/stepstone_small.parquet")
  rds <- testthat::test_path("fixtures/stepstone_small_targets.rds")
  skip_if_not(file.exists(pq))
  skip_if_not(file.exists(rds))
  ss  <- arrow::read_parquet(pq)
  tgt <- readRDS(rds)
  for (nm in names(tgt)) ss[[nm]] <- factor(ss[[nm]])
  # Force linear on a high-compression problem: will overflow → fallback to log
  withr::with_envvar(c(LBW_IEPPA_FORCE_PATH = "linear"), {
    r <- suppressWarnings(
      harvest(ss, tgt, method = "ieppa",
              max_iterations = 30L, accelerate = TRUE)
    )
  })
  res <- attr(r, "result")
  expect_true(is.finite(res$max_error),
              label = "fallback: finite max_error")
  expect_equal(sum(r$weights), nrow(ss), tolerance = 1e-6,
               label = "fallback: sum(weights) == n")
})
