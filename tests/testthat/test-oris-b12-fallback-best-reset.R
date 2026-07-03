context("oris/SRAA linear->log fallback best-iterate reset (CR-B12 / y2ks.12)")

## CR-B12: apply_linear_fallback reset best (BestIterTracker) and the SRAA error
## trackers, but NOT the SRAA best-iterate snapshot (sraa_has_best / lf_best).
## best.reset() clears best.has_best(), but the outer-stall revert (oris.cpp ~L1083)
## gates on sraa_has_best, which stayed true — so a post-fallback revert could
## restore the stale pre-fallback (degenerate linear-space) lf_best onto the log
## path, corrupting the returned solution. The fix clears sraa_has_best/lf_best in
## the fallback. Observable guard: after a forced linear->log fallback with SRAA,
## the returned solution stays valid (finite errors, Σw=n) — a restored degenerate
## iterate would break these.

test_that("forced linear->log fallback with accelerate returns a valid solution", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("withr")
  pq  <- testthat::test_path("fixtures/stepstone_small.parquet")
  rds <- testthat::test_path("fixtures/stepstone_small_targets.rds")
  skip_if_not(file.exists(pq))
  skip_if_not(file.exists(rds))
  ss  <- arrow::read_parquet(pq)
  tgt <- readRDS(rds)
  for (nm in names(tgt)) ss[[nm]] <- factor(ss[[nm]])
  n <- nrow(ss)
  # Force linear on a high-compression problem: linear overflows → apply_linear_fallback
  # → log path. A short budget then drives an outer-stall/BUDGET exit, exercising the
  # revert guard that reads sraa_has_best.
  withr::with_envvar(c(LBW_ORIS_FORCE_PATH = "linear"), {
    r <- suppressWarnings(
      harvest(ss, tgt, method = "oris", accelerate = TRUE, max_iterations = 40L)
    )
  })
  res <- attr(r, "result")

  expect_true(is.finite(res$max_error),
              label = "post-fallback: finite max_error (no degenerate revert)")
  expect_true(is.finite(res$best_error),
              label = "post-fallback: finite best_error")
  expect_gte(res$best_error, 0, label = "post-fallback: best_error non-negative")
  expect_equal(sum(r$weights), n, tolerance = 1e-6,
               label = "post-fallback: sum(weights) == n")
  # A stale degenerate iterate revert would leave a grossly mis-calibrated result;
  # the constrained error stays modest instead.
  expect_lt(res$max_error, 0.5,
            label = "post-fallback: constrained error stays sane")
})
