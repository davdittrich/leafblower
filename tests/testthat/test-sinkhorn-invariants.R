test_that("tsyw: sinkhorn Sigma_w=n holds under tight bounds", {
  skip_if_not_installed("leafblower")
  skip_if(!requireNamespace("arrow", quietly = TRUE), "arrow not available")
  # benchmarks/ lives at package root, two levels above tests/testthat/
  pkg_root   <- normalizePath(file.path(testthat::test_path(), "../.."))
  bench_data <- file.path(pkg_root, "benchmarks", "stepstone_bench_data.parquet")
  bench_tgt  <- file.path(pkg_root, "benchmarks", "stepstone_bench_targets.json")
  skip_if(!file.exists(bench_data), "stepstone bench data not available")
  skip_if(!file.exists(bench_tgt),  "stepstone bench targets not available")
  data <- arrow::read_parquet(bench_data)
  tgt_raw <- jsonlite::fromJSON(bench_tgt)
  target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })
  # min_weight=0.3 forces bisect_capacity_fast to clamp many cells to lower bound;
  # without at_lower[] tracking this causes Sigma_w != n at solver exit.
  suppressWarnings(
    w <- leafblower::harvest(data, target,
                             method = "sinkhorn",
                             max_weight = 5,
                             min_weight = 0.3,
                             max_iterations = 200L,
                             attach_weights = FALSE)
  )
  n <- nrow(data)
  expect_equal(sum(w), n, tolerance = 1e-6,
    label = paste0("sum(w)=", sum(w), " n=", n))
})
