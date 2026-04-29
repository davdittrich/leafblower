# tests/testthat/test-sraa-global.R
# T_sraa_global: SRAA-m global safeguard regression test (K=6 cross-margin overlap)

test_that("T_sraa_global: greenkhorn+SRAA max_err <= plain on K=6 cross-margin problem", {
  # K=6 with overlapping cross-margins reproduces the basin-escape failure.
  # Overlapping structure: gender in 3 margins (gender, gt, ga)
  #                        time   in 2 margins (time, gt, ta)
  #                        age    in 2 margins (age,  ga, ta)
  # RED with old per-step local safeguard (me_aa ≈ 7.9e-2 vs me_plain ≈ 3.7e-2, ratio ≈ 2.1).
  # GREEN after global best_err_seen safeguard + revert-to-best.
  set.seed(42); n <- 8000L
  df <- data.frame(
    gender = factor(sample(c("M", "F"), n, TRUE)),
    time   = factor(sample(1:3,         n, TRUE)),
    age    = factor(sample(1:4,         n, TRUE))
  )
  # Cross-margin interaction columns
  df$gt <- factor(paste0(as.character(df$gender), as.character(df$time)))
  df$ga <- factor(paste0(as.character(df$gender), as.character(df$age)))
  df$ta <- factor(paste0(as.character(df$time),   as.character(df$age)))

  gt_tbl <- table(df$gt) / n
  ga_tbl <- table(df$ga) / n
  ta_tbl <- table(df$ta) / n

  tgt <- list(
    gender = c(M = 0.48, F = 0.52),
    time   = setNames(c(0.4, 0.35, 0.25), as.character(1:3)),
    age    = setNames(c(0.3, 0.25, 0.25, 0.2), as.character(1:4)),
    gt     = { t <- setNames(as.numeric(gt_tbl) * c(0.95,1.02,0.98,1.03,0.97,1.05),
                             names(gt_tbl)); t / sum(t) },
    ga     = { t <- setNames(as.numeric(ga_tbl), names(ga_tbl)); t / sum(t) },
    ta     = { t <- setNames(as.numeric(ta_tbl), names(ta_tbl)); t / sum(t) }
  )
  tgt <- lapply(tgt, function(t) t / sum(t))

  r_aa    <- suppressWarnings(harvest(df, tgt, method = "greenkhorn", accelerate = TRUE,
                                       max_iterations = 150L, attach_weights = FALSE))
  r_plain <- suppressWarnings(harvest(df, tgt, method = "greenkhorn", accelerate = FALSE,
                                       max_iterations = 150L, attach_weights = FALSE))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label = sprintf(
      "SRAA K=6 cross-margin (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
})
