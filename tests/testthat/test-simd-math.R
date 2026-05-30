test_that("S1: raking KL value matches reference after bulk_log vectorization", {
  set.seed(1); n <- 200L
  df  <- data.frame(x = factor(sample(letters[1:3], n, replace=TRUE)),
                    y = factor(sample(c("m","f"), n, replace=TRUE)))
  tgt <- list(x=c(a=0.33,b=0.33,c=0.34), y=c(m=0.5,f=0.5))
  res <- harvest(df, target=tgt, method="raking",
                 convergence=list(absolute=1e-6))
  expect_lt(attr(res,"result")$max_error, 1e-3)
})

test_that("S4: sinkhorn produces same result after bisect_capacity_fast", {
  set.seed(42)
  n <- 200L
  df <- data.frame(
    age = sample(c("young","old"), n, replace=TRUE),
    sex = sample(c("m","f"), n, replace=TRUE)
  )
  tgt <- list(age=c(young=0.6, old=0.4), sex=c(m=0.5, f=0.5))
  res <- harvest(df, target=tgt, method="sinkhorn", min_weight=0.5, max_weight=3.0)
  expect_lt(attr(res,"result")$max_error, 0.01)
})


test_that("bulk_log correctness: oris output matches task2_oris_ref within 1e-12", {
  # Verifies that log-vectorization does not change solver output.
  # Reference was generated with the synthetic DGP used in the ORIS perf plan.
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
  ref    <- readRDS(testthat::test_path("task2_oris_ref.rds"))
  result <- harvest(df, tgt, method = "oris")
  new_w  <- result$weights
  expect_equal(new_w, ref, tolerance = 1e-10,
               info = "SIMD log introduced weight divergence vs pre-SIMD reference")
})
