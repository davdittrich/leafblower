# tests/testthat/test-sraa-global.R
# T_sraa_global: SRAA-m global safeguard regression test (K=4 overlapping margins)

test_that("T_sraa_global: greenkhorn+SRAA max_err <= plain on K=4 overlapping-margin problem", {
  # K=4 designed to reproduce the basin-escape failure seen at K=9 (stepstone).
  # Old local safeguard: AA overshoots max_err-optimal basin on multi-margin problems.
  # New global safeguard: AA stays in or returns to max_err-optimal basin.
  # RED with old per-step local safeguard (err_AA <= err_plain).
  # GREEN after global best_err_seen safeguard + revert-to-best.
  set.seed(5); n <- 3000L
  df <- data.frame(
    a = factor(sample(letters[1:4], n, TRUE)),
    b = factor(sample(LETTERS[1:3], n, TRUE)),
    c = factor(sample(c("x", "y"),  n, TRUE)),
    d = factor(sample(c("M", "F"),  n, TRUE))
  )
  tgt <- list(
    a = setNames(c(0.4, 0.3, 0.2, 0.1), letters[1:4]),
    b = setNames(c(0.5, 0.3, 0.2),      LETTERS[1:3]),
    c = c(x = 0.6, y = 0.4),
    d = c(M = 0.45, F = 0.55)
  )
  r_aa    <- suppressWarnings(harvest(df, tgt, method = "greenkhorn", accelerate = TRUE,
                                       max_iterations = 500L, attach_weights = FALSE))
  r_plain <- suppressWarnings(harvest(df, tgt, method = "greenkhorn", accelerate = FALSE,
                                       max_iterations = 500L, attach_weights = FALSE))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label = sprintf(
      "SRAA K=4 overlapping-margin (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
})

test_that("T_sraa_adaptive_K9: greenkhorn+AA K=9 max_err <= plain (stepstone)", {
  skip_if_not_installed("arrow"); skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  df$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
                function(t) { t <- unlist(t); t / sum(t) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  r_aa    <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=TRUE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  r_plain <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=FALSE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label=sprintf("K=9 SRAA (%.4e) must not exceed plain (%.4e)", me_aa, me_plain))
})

test_that("T_sraa_raking_K9: raking+AA K=9 max_err <= plain (stepstone)", {
  skip_if_not_installed("arrow"); skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  df$uuid <- NULL
  tgt <- lapply(jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
                function(t) { t <- unlist(t); t / sum(t) })
  for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  r_aa    <- suppressWarnings(harvest(df, tgt, method="raking", accelerate=TRUE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  r_plain <- suppressWarnings(harvest(df, tgt, method="raking", accelerate=FALSE,
                                       max_weight=5, min_weight=0, max_iterations=5000L,
                                       attach_weights=FALSE, verbose=0))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label=sprintf("K=9 raking+AA (%.4e) must not exceed plain (%.4e)", me_aa, me_plain))
})

test_that("T_sraa_outer_revert: greenkhorn+AA K=6 cross-margin quality <= plain", {
  # K=6 cross-margin problem that exhibits basin-escape. Combined fix must prevent it.
  set.seed(42); n <- 8000L
  df <- data.frame(gender=factor(sample(c("M","F"),n,TRUE)),
                   time=factor(sample(1:3,n,TRUE)), age=factor(sample(1:4,n,TRUE)))
  df$gt <- factor(paste0(df$gender,df$time))
  df$ga <- factor(paste0(df$gender,df$age))
  df$ta <- factor(paste0(df$time,df$age))
  gt_t <- table(df$gt)/n; ga_t <- table(df$ga)/n; ta_t <- table(df$ta)/n
  tgt <- list(
    gender=c(M=0.48,F=0.52), time=setNames(c(0.4,0.35,0.25),1:3),
    age=setNames(c(0.3,0.25,0.25,0.2),1:4),
    gt={t<-setNames(as.numeric(gt_t)*c(0.95,1.02,0.98,1.03,0.97,1.05),names(gt_t));t/sum(t)},
    ga={t<-setNames(as.numeric(ga_t),names(ga_t));t/sum(t)},
    ta={t<-setNames(as.numeric(ta_t),names(ta_t));t/sum(t)})
  tgt <- lapply(tgt, function(t) t/sum(t))
  r_aa    <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=TRUE,
                                       max_iterations=200L,attach_weights=FALSE))
  r_plain <- suppressWarnings(harvest(df,tgt,method="greenkhorn",accelerate=FALSE,
                                       max_iterations=200L,attach_weights=FALSE))
  me_aa    <- attr(r_aa,   "result")$max_error
  me_plain <- attr(r_plain,"result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label=sprintf("K=6 combined fix: AA (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
})
