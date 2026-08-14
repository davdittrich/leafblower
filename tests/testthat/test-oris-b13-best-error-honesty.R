## CR-B13: on the ORIS accelerate (SRAA) path, res$best_error was copied from
## best.best_metric, tracked on the PRE-clamp accepted cell masses. On a
## tight/infeasible exit the returned bounds-clamped weights miss margins (large
## error) yet best_error read ~ -1.11e-16 ("converged"). harvest.R:728,740 surfaces
## best_error verbatim in the status==4L BUDGET warning, LABELLED with
## convergence_used$metric. The fix recomputes best_error on the finalized clamped
## weights in the CONVERGENCE-METRIC family (select_metric(metric, cm)), reusing the
## same CellMetrics mxcl.6 aggregates for max_error (no second aggregation). It is
## honest (magnitude of the returned solution) and correctly labelled (its units
## match convergence_used$metric and metric_prev_check), not the pre-clamp ~0 lie.

test_that("accelerate best_error is honest (metric-family, not the pre-clamp lie)", {
  set.seed(99L)
  n <- 1000L
  df <- data.frame(
    x = factor(sample(c("H", "L"), n, replace = TRUE, prob = c(0.7, 0.3))),
    y = factor(sample(c("P", "Q"), n, replace = TRUE, prob = c(0.3, 0.7)))
  )
  tgt <- list(x = c(H = 0.4, L = 0.6), y = c(P = 0.7, Q = 0.3))
  # Conflicting margins + tight bounds → constrained optimum with residual error.
  r <- suppressWarnings(
    harvest(df, tgt, method = "oris", accelerate = TRUE,
            max_weight = 1.5, min_weight = 0.5, max_iterations = 80L)
  )
  res <- attr(r, "result")

  expect_true(is.finite(res$best_error), label = "best_error is finite")
  # Pre-fix this read ~ -1.11e-16 (a nonsensical negative "converged" value).
  expect_gte(res$best_error, 0, label = "best_error non-negative (not the pre-clamp lie)")
  # Fixture genuinely misses margins → the returned solution has real error.
  expect_gt(res$max_error, 1e-3, label = "fixture reaches a constrained optimum")
  # Honest: the post-clamp constrained metric of the RETURNED solution far exceeds
  # the pre-clamp mid-loop value (metric_prev_check) that fed the false ~0 reading.
  expect_gt(res$best_error, 1e-2,
            label = "best_error reflects the constrained solution, not ~0")
  expect_gt(res$best_error, res$metric_prev_check,
            label = "post-clamp best_error > pre-clamp metric_prev_check")
})

test_that("accelerate BUDGET-exit best_error is honest at the harvest warning surface", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("withr")
  pq  <- testthat::test_path("fixtures/stepstone_small.parquet")
  rds <- testthat::test_path("fixtures/stepstone_small_targets.rds")
  skip_if_not(file.exists(pq))
  skip_if_not(file.exists(rds))
  ss  <- arrow::read_parquet(pq)
  tgt <- readRDS(rds)
  for (nm in names(tgt)) ss[[nm]] <- factor(ss[[nm]])
  # Forced-linear on a high-compression problem overflows → log fallback, then
  # exhausts a short budget → status==4L (BUDGET), the harvest.R:727 warning surface.
  withr::with_envvar(c(LBW_ORIS_FORCE_PATH = "linear"), {
    r <- suppressWarnings(
      harvest(ss, tgt, method = "oris", accelerate = TRUE, max_iterations = 40L)
    )
  })
  res <- attr(r, "result")
  expect_equal(res$status, 4L, label = "status is BUDGET (4L)")
  expect_true(is.finite(res$best_error), label = "best_error finite")
  expect_gte(res$best_error, 0, label = "best_error non-negative")
  # Honest constrained metric of the returned solution, not the pre-clamp ~0 lie.
  expect_gt(res$best_error, 1e-5, label = "best_error is the honest returned-solution error")
})

test_that("oris BUDGET exit emits the generic budget message, never a geometric-rate projection (hyfk / CR-B13b)", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("withr")
  pq  <- testthat::test_path("fixtures/stepstone_small.parquet")
  rds <- testthat::test_path("fixtures/stepstone_small_targets.rds")
  skip_if_not(file.exists(pq))
  skip_if_not(file.exists(rds))
  ss  <- arrow::read_parquet(pq)
  tgt <- readRDS(rds)
  for (nm in names(tgt)) ss[[nm]] <- factor(ss[[nm]])
  msgs <- character(0)
  withr::with_envvar(c(LBW_ORIS_FORCE_PATH = "linear"), {
    r <- withCallingHandlers(
      harvest(ss, tgt, method = "oris", accelerate = TRUE, max_iterations = 40L,
              attach_weights = FALSE),
      warning = function(cond) { msgs <<- c(msgs, conditionMessage(cond)); invokeRestart("muffleWarning") })
  })
  res <- attr(r, "result")
  # Forced-linear on this high-compression fixture reliably reaches BUDGET (the
  # sibling test above asserts status==4L unconditionally on the same inputs), so
  # assert it directly rather than skipping — a silent skip would hide coverage.
  expect_equal(res$status, 4L, label = "forced-linear oris reaches BUDGET (4L)")
  budget_msg <- grep("budget exhausted", msgs, value = TRUE)
  expect_length(budget_msg, 1L)
  # oris is box-constrained water-fill: under active clamps it converges
  # piecewise-linearly / slow-rate (O(t^-1/2)), so the geometric-rate projection
  # ("Asymptotic rate ... total iterations needed") is model-invalid and must
  # NEVER be emitted for oris. The generic "Increase max_iterations" is honest.
  expect_false(grepl("Asymptotic rate", budget_msg),
               label = "oris BUDGET must not emit a geometric-rate projection")
  expect_true(grepl("Increase max_iterations", budget_msg),
              label = "oris BUDGET emits the generic budget message")
})
