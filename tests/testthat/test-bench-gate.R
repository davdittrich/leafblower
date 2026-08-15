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
  expect_lte(r$max_w, r$max_weight + 1e-10, label = "max_w: per-observation bound honoured")
  expect_gte(r$min_w, -1e-10, label = "min_w: no negative weight produced")
  expect_lte(r$max_error, 1e-3, label = "max_error: survey-adequate accuracy floor")
  expect_true(is.finite(r$wall_s) && r$wall_s > 0, label = "wall_s: finite and positive")
  # Paired headline claim (D-06 read with D-10, developer-selected
  # 2026-08-15): speed is never asserted without accuracy/bound/n_eff
  # alongside it. Both ceilings below apply to
  # leafblower_oris_soft/medium_100k_5margins.
  #
  # wall_s <= 0.5s: derived as half of the slowest doc-named competitor's
  # measured time on this class (ReGenesees_e_calibrate, 0.9051s, from
  # benchmarks/results/oris_soft_vs_competitors.csv), giving ~11.7x
  # headroom over the measured 0.0427s while staying meaningfully faster
  # than every competitor even under adverse host conditions. Machine:
  # AMD Ryzen 9 9950X3D 16-Core Processor, R 4.6.1 (2026-06-24),
  # single-thread BLAS/LAPACK (see
  # benchmarks/results/oris_soft_vs_competitors_env.txt). Wall time is the
  # most machine-sensitive column in this table, so the ceiling is set with
  # deliberate headroom rather than at the measured value -- a gate pinned
  # to the measurement would fail on a busier host and get disabled.
  expect_lte(r$wall_s, 0.5,
             label = "wall_s: <=0.5s ceiling on leafblower_oris_soft/medium_100k_5margins")
  # n_eff >= 60000: measured n_eff=67489.4 on this row (same CSV), ~11%
  # headroom below the measurement -- the same margin the pre-existing
  # max_error<=1e-3 floor keeps below its own measured 3.35e-05. A floor
  # clearly cleared today, sized to catch a real accuracy/deff regression
  # rather than host noise.
  expect_gte(r$n_eff, 60000,
             label = "n_eff: >=60000 floor on leafblower_oris_soft/medium_100k_5margins")
  cat(sprintf("honest gate: wall_s=%.4f max_error=%.3e max_w=%.4f n_eff=%.1f\n",
              r$wall_s, r$max_error, r$max_w, r$n_eff))
})

# kk1204-shape regression floor (D-01/D-11, re-gated 2026-08-15): this fixture shares
# n=500000/K=20/seed=1204/max_weight=3 with docs/performance.md's documented
# known_limit_k20_uniform class, but uses a literal uniform 1/5 target per margin, NOT
# the skewed target (0.3/0.175x4) that reproduces the near-infeasible degenerate case
# documented in docs/performance.md's "Known limit" section and
# docs/investigations/2026-04-23-kk1204-convergence.md (commit 3effd3a). Under the
# uniform target this block converges trivially (see pre-edit observation below), so it
# is NOT a test of the degenerate known limit -- it is a shape-matched (n/K) regression
# floor confirming oris still converges fast and accurately at this size. D-01 dropped
# this class as the basis of any headline claim; leafblower-ylsy (closed) is cited, not
# reopened, per D-03.
#
# Pre-edit observation (2026-08-15, before this block's guard/label rewrite, CI unset,
# NOT_CRAN=true): status=0, iterations=10, best_error=-7.376e-14, elapsed=1.6s -- PASS
# on all three assertions below. The existing best_error<=1e-3 and elapsed<=30
# thresholds are reconfirmed against this observation with large headroom (both cleared
# by several orders of magnitude); this is NOT the 03-02 known_limit_k20_uniform
# measurement (max_error 5.229e-03 on the skewed target) -- that is a different fixture
# and its numbers are not comparable to this one's.
test_that("kk1204-shape regression floor: n=500k K=20 uniform target converges in <30s with best_error<1e-3", {
  skip_on_cran()
  skip_if(Sys.getenv("LBW_BENCH_GATE") == "")
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
