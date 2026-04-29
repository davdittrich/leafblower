# tests/testthat/test-sraa-global.R
# T_sraa_global: SRAA-m global safeguard regression test

# ── T1: K=4 independent margins (original test, must stay GREEN) ────────────
test_that("T_sraa_global: greenkhorn+SRAA max_err <= plain on K=4 overlapping-margin problem", {
  # K=4 chosen to reproduce the basin-escape failure seen at K=9 (stepstone).
  # Old local safeguard: AA overshoots max_err-optimal basin on multi-margin problems.
  # New global safeguard: AA stays in or returns to max_err-optimal basin.
  set.seed(5); n <- 3000L
  df <- data.frame(
    a = factor(sample(letters[1:4], n, TRUE)),
    b = factor(sample(LETTERS[1:3], n, TRUE)),
    c = factor(sample(c("x","y"),   n, TRUE)),
    d = factor(sample(c("M","F"),   n, TRUE))
  )
  tgt <- list(
    a = setNames(c(0.4, 0.3, 0.2, 0.1), letters[1:4]),
    b = setNames(c(0.5, 0.3, 0.2),      LETTERS[1:3]),
    c = c(x=0.6, y=0.4),
    d = c(M=0.45, F=0.55)
  )
  r_aa    <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=TRUE,
                                       max_iterations=500L, attach_weights=FALSE))
  r_plain <- suppressWarnings(harvest(df, tgt, method="greenkhorn", accelerate=FALSE,
                                       max_iterations=500L, attach_weights=FALSE))
  me_aa    <- attr(r_aa,    "result")$max_error
  me_plain <- attr(r_plain, "result")$max_error
  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label=sprintf("SRAA K=4 (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
})

# ── T2: K=6 cross-margin interactions — RED gate for global safeguard ────────
# Basin-escape manifests when shared factors appear across multiple margins:
#   gender in 3 margins (gender, gt, ga)
#   time   in 2 margins (time, gt, ta)
#   age    in 2 margins (age, ga, ta)
# max_iterations=150 keeps convergence at a visible plateau (1e-2 range), not
# machine epsilon, so the ratio (>2x) is unmistakable.
#
# RED: greenkhorn+SRAA throws an error and reports max_error=7.91e-02
#       while plain converges to max_error=3.72e-02 — ratio ≈ 2.13.
# Stays RED until a global safeguard (tracking best-ever weights) is implemented.
test_that("T_sraa_global: greenkhorn+SRAA max_err <= plain on K=6 cross-margin problem", {
  set.seed(42); n <- 8000L
  df <- data.frame(
    gender = factor(sample(c("M", "F"), n, TRUE)),
    time   = factor(sample(1:3,         n, TRUE)),
    age    = factor(sample(1:4,         n, TRUE))
  )
  # Cross-margin interaction columns (no sep=":" — avoids C-bridge NULL issue)
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

  # greenkhorn+SRAA throws an error on divergence; catch and parse max_error.
  # The error message format is: "leafblower: greenkhorn: N iters, max_error=X.XXe-XX"
  me_aa <- tryCatch({
    r <- suppressWarnings(
      harvest(df, tgt, method = "greenkhorn", accelerate = TRUE,
              max_iterations = 150L)
    )
    attr(r, "result")$max_error
  }, error = function(e) {
    m <- regmatches(conditionMessage(e),
                    regexpr("max_error=([0-9eE.+\\-]+)", conditionMessage(e)))
    if (length(m) == 0L)
      stop("unexpected error in greenkhorn+SRAA: ", conditionMessage(e))
    as.numeric(sub("max_error=", "", m))
  })

  r_plain <- suppressWarnings(
    harvest(df, tgt, method = "greenkhorn", accelerate = FALSE,
            max_iterations = 150L)
  )
  me_plain <- attr(r_plain, "result")$max_error

  expect_lte(me_aa, me_plain * 1.001 + 1e-10,
    label = sprintf(
      "SRAA K=6 cross-margin (%.2e) must not exceed plain (%.2e)", me_aa, me_plain))
})
