test_that("method='rake' emits warning about L-BFGS-B", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a","b"), 200, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))
  expect_warning(harvest(df, tgt, method="rake", convergence = list(absolute = 1e-6)), regexp = "L-BFGS-B")
})

test_that("method='nr' emits warning about L-BFGS-B", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a","b"), 200, replace=TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5))
  expect_warning(harvest(df, tgt, method="nr", convergence = list(absolute = 1e-6)), regexp = "L-BFGS-B")
})

test_that("default routing selects iEPPA for large complexity", {
  set.seed(1)
  n   <- 200000L
  df  <- data.frame(x = factor(sample(c("a","b","c"), n, replace=TRUE)))
  tgt <- list(x = c(a=0.33, b=0.34, c=0.33))
  # iEPPA is the default; no method arg needed
  result <- harvest(df, tgt, convergence = list(absolute = 1e-6))
  expect_identical(attr(result, "algorithm"), "ieppa")
})

test_that("convergence$absolute is forwarded to solver", {
  set.seed(42)
  n <- 10000L
  df <- data.frame(
    age = factor(sample(c("Y","M","O"), n, replace=TRUE,
                        prob=c(0.60, 0.30, 0.10))),
    sex = factor(sample(c("M","F"),     n, replace=TRUE,
                        prob=c(0.70, 0.30)))
  )
  tgt <- list(
    age = c(Y=0.33, M=0.40, O=0.27),
    sex = c(M=0.49, F=0.51)
  )
  # 2 iterations with default tol (1e-6): competing margins cannot converge -> warning
  expect_warning(
    harvest(df, tgt, method="ieppa", max_iterations=2, convergence = list(absolute = 1e-6)),
    regexp="did not converge"
  )
  # 2 iterations with loose tol (0.3): error after 2 iters < 0.3 -> no warning
  # Before fix: tol_abs ignored, tol=1e-6 used, warning fires -> test fails
  # After fix:  tol_abs=0.3 forwarded, error < 0.3 accepted -> no warning
  expect_no_warning(
    harvest(df, tgt, method="ieppa", max_iterations=2,
            convergence=list(absolute=0.3))
  )
})

test_that("normalize_start_weights rejects all-zero vector", {
  df  <- data.frame(x = factor(c("a","b","a")))
  tgt <- list(x = c(a=0.5, b=0.5))
  expect_error(
    harvest(df, tgt, start_weights = c(0, 0, 0), convergence = list(absolute = 1e-6)),
    regexp = "start_weights must sum to a positive value"
  )
})

test_that("parse_target stops on unnamed 3-column data frame", {
  tgt_df <- data.frame(v = c("x","x"), l = c("a","b"), p = c(0.5, 0.5))
  df     <- data.frame(x = factor(c("a","b")))
  expect_error(
    harvest(df, tgt_df, convergence = list(absolute = 1e-6)),
    regexp = "no 'variable'/'level'/'proportion' names"
  )
})

test_that("infeasible problem (empty cell with positive target) stops with infeasible error", {
  # All 3 observations are in group "a"; no observation is in group "b".
  # Target for "b" = 0.4 > 0 → empty cell with positive target → infeasible.
  df  <- data.frame(x = factor(c("a", "a", "a"), levels = c("a", "b")))
  tgt <- list(x = c(a = 0.6, b = 0.4))
  expect_error(
    harvest(df, tgt, method = "ieppa", convergence = list(absolute = 1e-6)),
    regexp = "infeasible problem"
  )
})

test_that("K > 64 rejected with informative error", {
  # 65 margin columns — should fail validation
  n <- 100
  data <- as.data.frame(matrix("0", nrow = n, ncol = 65))
  names(data) <- paste0("v", 1:65)
  for (i in seq_len(65)) data[[i]] <- factor(data[[i]])
  # Targets: each column has 1 category, target = 1
  tgt <- lapply(seq_len(65), function(k) setNames(1.0, "0"))
  names(tgt) <- names(data)
  expect_error(harvest(data, tgt, method = "ieppa", convergence = list(absolute = 1e-6)),
               regexp = "K.*64|too many margin", ignore.case = TRUE)
})

test_that("cat_counts <= 0 rejected", {
  # Empty target list for a column — degenerate
  n <- 100
  data <- data.frame(a = sample(c("x", "y"), n, replace = TRUE))
  tgt <- list(a = setNames(numeric(0), character(0)))
  expect_error(harvest(data, tgt, method = "ieppa", convergence = list(absolute = 1e-6)),
               regexp = "cat_counts|empty target|no categories", ignore.case = TRUE)
})

test_that("zero-sum input weights rejected", {
  n <- 100
  data <- data.frame(a = sample(c("x", "y"), n, replace = TRUE))
  tgt <- list(a = c(x = 0.5, y = 0.5))
  sw <- rep(0.0, n)
  expect_error(harvest(data, tgt, method = "ieppa", start_weights = sw,
                       convergence = list(absolute = 1e-6)),
               regexp = "start_weights.*positive|sum.*zero", ignore.case = TRUE)
})
