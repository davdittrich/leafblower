test_that("bulk_log correctness: ieppa output matches task2_ieppa_ref within 1e-12", {
  # Verifies that log-vectorization does not change solver output.
  # Reference was generated with the synthetic DGP used in the iEPPA perf plan.
  set.seed(42)
  n   <- 10000L
  df  <- data.frame(
    age = factor(sample(c("18-34", "35-54", "55+"), n, replace = TRUE)),
    sex = factor(sample(c("M", "F"), n, replace = TRUE))
  )
  tgt <- list(
    age = c("18-34" = 0.30, "35-54" = 0.45, "55+" = 0.25),
    sex = c(M = 0.50, F = 0.50)
  )
  ref    <- readRDS(testthat::test_path("task2_ieppa_ref.rds"))
  result <- harvest(df, tgt, method = "ieppa")
  new_w  <- result$weights
  expect_equal(new_w, ref, tolerance = 1e-10,
               info = "SIMD log introduced weight divergence vs pre-SIMD reference")
})
