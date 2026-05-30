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
