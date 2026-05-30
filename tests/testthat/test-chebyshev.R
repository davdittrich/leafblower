# tests/testthat/test-chebyshev.R
# RED tests for chebyshev Mehrotra warm-start rewrite + greg quality warning

test_that("T_greg_warn: greg warns when max_err > 5% on K=5 tight-bounds problem", {
  set.seed(99); n <- 3000L
  df <- data.frame(
    a = factor(sample(letters[1:4], n, TRUE)),
    b = factor(sample(LETTERS[1:3], n, TRUE)),
    c = factor(sample(c("x","y"),   n, TRUE)),
    d = factor(sample(c("M","F"),   n, TRUE)),
    e = factor(sample(c("Y","O"),   n, TRUE))
  )
  tgt <- list(
    a = setNames(c(0.4,0.3,0.2,0.1), letters[1:4]),
    b = setNames(c(0.5,0.3,0.2),     LETTERS[1:3]),
    c = c(x=0.6,y=0.4), d = c(M=0.45,F=0.55), e = c(Y=0.55,O=0.45)
  )
  expect_warning(
    harvest(df, tgt, method="greg", max_weight=1.8, min_weight=0,
            max_iterations=50L, attach_weights=FALSE),
    regexp="greg.*unreliable|greg.*max_err",
    label="greg must warn on K=5 tight-bounds (max_err expected >5%)")
})

test_that("T_cheby_warm: chebyshev K=3 converges with status=0 and max_err <= raking", {
  set.seed(7); n <- 5000L
  df <- data.frame(
    a = factor(sample(letters[1:3], n, TRUE)),
    b = factor(sample(LETTERS[1:4], n, TRUE)),
    c = factor(sample(c("M","F"),   n, TRUE))
  )
  tgt <- list(
    a = setNames(c(0.4,0.35,0.25), letters[1:3]),
    b = setNames(c(0.3,0.3,0.2,0.2), LETTERS[1:4]),
    c = c(M=0.48, F=0.52)
  )
  r_cheby  <- suppressWarnings(harvest(df, tgt, method="chebyshev",
                                        max_iterations=500L, attach_weights=FALSE))
  r_raking <- suppressWarnings(harvest(df, tgt, method="raking",
                                        max_iterations=500L, attach_weights=FALSE))
  me_c <- attr(r_cheby,  "result")$max_error
  me_r <- attr(r_raking, "result")$max_error
  st_c <- attr(r_cheby,  "result")$status
  expect_equal(st_c, 0L,
    label="chebyshev must converge (status=0) on K=3")
  expect_true(is.finite(me_c),
    label="chebyshev max_error must be finite")
  expect_lte(me_c, me_r * 1.001 + 1e-10,
    label=sprintf("chebyshev (%.4e) must not exceed raking (%.4e)", me_c, me_r))
  expect_lt(me_c, 1e-3,
    label="chebyshev must converge to <1e-3 on K=3")
})

test_that("T_cheby_warm_fallback: chebyshev returns finite result even with marginal warm-start", {
  # With very few iterations, oris warm-start is weak but chebyshev should
  # still return a finite (if not converged) result without crashing.
  set.seed(7); n <- 1000L
  df  <- data.frame(x=factor(sample(letters[1:3],n,TRUE)),
                    y=factor(sample(c("M","F"),n,TRUE)))
  tgt <- list(x=c(a=0.3,b=0.4,c=0.3), y=c(M=0.5,F=0.5))
  r <- suppressWarnings(harvest(df, tgt, method="chebyshev",
                                 max_iterations=200L, attach_weights=FALSE))
  expect_true(is.finite(attr(r, "result")$max_error),
    label="chebyshev must return finite max_error (not NaN from bad warm-start)")
})

test_that("T_cheby_K4: chebyshev converges on K=4 overlapping-margin problem (ν fix)", {
  set.seed(7L); n <- 2000L
  df <- data.frame(
    a = factor(sample(3L, n, replace = TRUE)),
    b = factor(sample(4L, n, replace = TRUE)),
    c = factor(sample(3L, n, replace = TRUE)),
    d = factor(sample(2L, n, replace = TRUE))
  )
  tgt <- list(
    a = setNames(rep(1/3, 3L), as.character(1:3)),
    b = setNames(rep(1/4, 4L), as.character(1:4)),
    c = setNames(rep(1/3, 3L), as.character(1:3)),
    d = setNames(rep(1/2, 2L), as.character(1:2))
  )
  r <- suppressWarnings(
    harvest(df, tgt, method = "chebyshev", max_iterations = 300L,
            attach_weights = FALSE, verbose = 0L)
  )
  res <- attr(r, "result")
  expect_equal(res$status, 0L,
    label = sprintf("K=4 chebyshev: status=%d max_err=%.2e", res$status, res$max_error))
  expect_lt(res$max_error, 1e-3,
    label = sprintf("K=4 chebyshev: max_err=%.2e", res$max_error))
})

test_that("T_cheby_K9: chebyshev K=9 stepstone max_err <= greenkhorn", {
  skip_if_not_installed("arrow")
  skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  df$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
                function(t) { t <- unlist(t); t / sum(t) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  r_c <- suppressWarnings(harvest(df, tgt, method="chebyshev",
                                   max_weight=5, min_weight=0, max_iterations=5000L,
                                   attach_weights=FALSE, verbose=0))
  r_g <- suppressWarnings(harvest(df, tgt, method="greenkhorn",
                                   max_weight=5, min_weight=0, max_iterations=5000L,
                                   attach_weights=FALSE, verbose=0))
  me_c <- attr(r_c, "result")$max_error
  me_g <- attr(r_g, "result")$max_error
  expect_equal(attr(r_c, "result")$status, 0L,
    label="chebyshev must converge on K=9")
  expect_lte(me_c, me_g * 1.001 + 1e-10,
    label=sprintf("chebyshev (%.4e) must not exceed greenkhorn (%.4e)", me_c, me_g))
})
