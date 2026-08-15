test_that("stepstone-fulldata AB config meets merge floor (errRp) + Pearson agreement", {
  skip_on_cran()
  # Only run when explicitly triggered (not in normal devtools::test() runs).
  skip_if(Sys.getenv("LBW_BENCH_GATE") == "")
  rpt_path <- "benchmarks/stepstone_fulldata_homotopy_report.rds"
  skip_if(!file.exists(rpt_path))
  report   <- readRDS(rpt_path)
  ab_row   <- report[report$config == "leafblower_AB", ]
  # Gate uses internal errRp (comparable to lb_max_error=0.00222 in reference commit 9a97cc8).
  # tapply max_err and errRp are different metrics — gate uses errRp column.
  errRp_col <- if ("errRp" %in% names(report)) "errRp" else "max_err"
  expect_lte(ab_row[[errRp_col]], 1.60e-3)
  pearson <- attr(report, "pearson_r")
  if (!is.null(pearson) && !is.na(pearson[["AB"]])) {
    expect_gte(pearson[["AB"]], 0.99)
  }
})

test_that("stepstone-fulldata rate slope <= -0.75 (beats O(t^-0.5) baseline)", {
  skip_on_cran()
  skip_if(Sys.getenv("LBW_BENCH_GATE") == "")
  rds_path <- "benchmarks/final_rate_slope.rds"
  skip_if(!file.exists(rds_path))
  res <- readRDS(rds_path)
  expect_lte(res$slope, -0.75)
})

test_that("honest gate: leafblower_oris_soft medium_100k_5margins bound/accuracy", {
  skip_on_cran()
  skip_if(Sys.getenv("LBW_BENCH_GATE") == "")
  # Bare-relative "benchmarks/..." (the sibling tests' idiom above) breaks
  # under a filtered testthat::test_dir(filter="bench-gate") run: only
  # test-algo-selection.R's own top-level setwd() resets cwd to the package
  # root, and a filter excluding that file never executes it, leaving cwd at
  # tests/testthat/ for the whole run. testthat::test_path() is cwd-agnostic
  # by design (test-sinkhorn-invariants.R already anchors on it the same
  # way) and resolves correctly under both devtools::test() and a filtered
  # test_dir() call, so it is used here instead of a bare relative string.
  pkg_root <- normalizePath(file.path(testthat::test_path(), "../.."))
  csv_path <- file.path(pkg_root, "benchmarks", "results", "oris_soft_vs_competitors.csv")
  skip_if(!file.exists(csv_path))
  d <- read.csv(csv_path)
  r <- d[d$arm == "leafblower_oris_soft" & d$input_class == "medium_100k_5margins", ]
  expect_equal(nrow(r), 1L, label = "exactly one leafblower_oris_soft/medium_100k_5margins row")
  # Wall-time ceiling deliberately NOT asserted here: the published headline
  # number and this gate's threshold are a one-way decision (once in the
  # README it cannot be quietly lowered) set in plan 03 after a human sees
  # the measured values (D-06 read with D-10). Wired in plan 03, not
  # forgotten here.
  expect_lte(r$max_w, r$max_weight + 1e-10, label = "max_w: per-observation bound honoured")
  expect_gte(r$min_w, -1e-10, label = "min_w: no negative weight produced")
  expect_lte(r$max_error, 1e-3, label = "max_error: survey-adequate accuracy floor")
  expect_true(is.finite(r$wall_s) && r$wall_s > 0, label = "wall_s: finite and positive")
  cat(sprintf("honest gate: wall_s=%.4f max_error=%.3e max_w=%.4f n_eff=%.1f\n",
              r$wall_s, r$max_error, r$max_w, r$n_eff))
})

test_that("kk1204 gate: n=500k K=20 converges in <30s with best_error<1e-3", {
  skip_on_cran()
  skip_if(Sys.getenv("CI") != "")
  set.seed(1204)
  n <- 500000; K <- 20; cats <- 5
  data_kk <- as.data.frame(
    replicate(K, factor(sample(seq_len(cats), n, replace = TRUE)),
              simplify = FALSE))
  names(data_kk) <- paste0("m", seq_len(K))
  target_kk <- stats::setNames(
    lapply(seq_len(K), function(k)
      setNames(rep(1/cats, cats), as.character(seq_len(cats)))),
    names(data_kk))
  t0 <- proc.time()["elapsed"]
  w <- leafblower::harvest(
    data_kk, target_kk, max_weight = 3,
    method = "oris",
    max_iterations = 500,
    attach_weights = FALSE)
  elapsed <- proc.time()["elapsed"] - t0
  r <- attr(w, "result")
  cat(sprintf("kk1204 gate: status=%d iters=%d best_error=%.3e time=%.1fs\n",
              r$status, r$iterations, r$best_error, elapsed))
  expect_equal(r$status, 0L, label = "status: must converge with default max_err+improvement criterion")
  expect_lte(r$best_error, 1e-3, label = "best_error: quality gate below 1e-3")
  expect_lte(elapsed, 30, label = "elapsed: speed gate <30s on n=500k K=20")
})
