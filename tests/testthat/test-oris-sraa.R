context("oris SRAA-m acceleration")

# Shared multi-margin fixture that uses the linear path.
# stepstone_small: n=10000, K=9, M_cell~5980, compression=1.7x → path=linear.
# make_tight_k5 had K=1 → log path → SRAA inactive. Use stepstone_small instead.
load_stepstone <- function() {
  skip_if_not_installed("arrow")
  pq  <- testthat::test_path("fixtures/stepstone_small.parquet")
  rds <- testthat::test_path("fixtures/stepstone_small_targets.rds")
  skip_if_not(file.exists(pq))
  skip_if_not(file.exists(rds))
  ss  <- arrow::read_parquet(pq)
  tgt <- readRDS(rds)
  for (nm in names(tgt)) ss[[nm]] <- factor(ss[[nm]])
  list(df = ss, tgt = tgt)
}

## Test 1 — SRAA activates and accepts AA steps ----------------------------

test_that("oris accelerate=TRUE activates SRAA and accepts AA steps (aa_count > 0)", {
  fx <- load_stepstone()
  r  <- suppressWarnings(
    harvest(fx$df, fx$tgt, method = "oris", max_iterations = 200L, accelerate = TRUE)
  )
  cnt <- attr(r, "result")$aa_accepted_count
  expect_gte(cnt, 1L,
             label = sprintf("aa_accepted_count = %d, expected >= 1", cnt))
})

## Test 2 — Regression: accel produces valid weights at equal budget --------
# oris minimizes marginal KL, not max_err. SRAA can produce higher max_err
# than plain when it converges to a different KL-optimal point. The
# regression guard checks validity (finite max_error, correct weight sum)
# rather than max_err parity.

test_that("oris accelerate=TRUE produces valid weights at equal budget (regression guard)", {
  fx     <- load_stepstone()
  r_accel <- suppressWarnings(
    harvest(fx$df, fx$tgt, method = "oris", max_iterations = 100L, accelerate = TRUE)
  )
  res <- attr(r_accel, "result")
  expect_true(is.finite(res$max_error),
              label = "accel max_error is finite")
  expect_true(res$status %in% c(0L, 1L, 4L, 5L),
              label = sprintf("accel status=%d is valid", res$status))
  expect_equal(sum(r_accel$weights), nrow(fx$df), tolerance = 1e-6,
               label = "sum(weights) == n")
})

## Test 3 — oris_soft runs with accelerate=TRUE ---------------------------

test_that("oris_soft accelerate=TRUE runs without error on stepstone_small", {
  fx <- load_stepstone()
  r  <- suppressWarnings(
    harvest(fx$df, fx$tgt, method = "oris_soft",
            max_weight = 5.0, min_weight = 0.0,
            max_iterations = 100L, accelerate = TRUE)
  )
  res <- attr(r, "result")
  expect_true(res$status %in% c(0L, 1L, 4L, 5L),
              label = sprintf("oris_soft + accel: status=%d", res$status))
  expect_true(is.finite(res$max_error),
              label = "oris_soft + accel: finite max_error")
})

## Test 4 — SRAA active on log path with accelerate=TRUE --------------------

test_that("oris accelerate=TRUE activates SRAA on log path at verbose=1", {
  # LL3: sraa_active_lvl = st.accelerate (path-agnostic). This K=2 case uses
  # the log path (compression=50x); SRAA now runs and logs [sraa] lines.
  set.seed(1L); n <- 200L
  df  <- data.frame(x = factor(sample(c("a","b"), n, TRUE)),
                    y = factor(sample(c("p","q"), n, TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5), y = c(p=0.4, q=0.6))
  r <- suppressWarnings(
    harvest(df, tgt, method = "oris", scheduler = "greedy",
            accelerate = TRUE, max_iterations = 20L, verbose = 1L)
  )
  # SRAA logs go via Rprintf (stdout), not R messages — check aa_accepted_count.
  expect_gt(attr(r, "result")$aa_accepted_count, 0L,
            label = "aa_accepted_count > 0 on log path")
})

## Test 5 — Output correlation with plain ----------------------------------

test_that("oris accelerate=TRUE produces highly correlated weights with plain", {
  fx  <- load_stepstone()
  r_p <- suppressWarnings(
    harvest(fx$df, fx$tgt, method = "oris", max_iterations = 200L, accelerate = FALSE)
  )
  r_a <- suppressWarnings(
    harvest(fx$df, fx$tgt, method = "oris", max_iterations = 200L, accelerate = TRUE)
  )
  expect_gt(cor(r_p$weights, r_a$weights), 0.95,
            label = "cor(plain, accel) > 0.95")
})

## Test 6 — aa_accepted_count field accessible from R ---------------------

test_that("aa_accepted_count field is accessible via attr(r, 'result')$aa_accepted_count", {
  fx  <- load_stepstone()
  r_a <- suppressWarnings(
    harvest(fx$df, fx$tgt, method = "oris", max_iterations = 200L, accelerate = TRUE)
  )
  r_p <- suppressWarnings(
    harvest(fx$df, fx$tgt, method = "oris", max_iterations = 200L, accelerate = FALSE)
  )
  cnt_a <- attr(r_a, "result")$aa_accepted_count
  cnt_p <- attr(r_p, "result")$aa_accepted_count
  expect_true(is.integer(cnt_a) || is.numeric(cnt_a),
              label = "aa_accepted_count is numeric")
  expect_equal(cnt_p, 0L,
               label = "aa_accepted_count = 0 when accelerate=FALSE")
})

## Test 7 — Budget regression: accel never catastrophically worse ---------

test_that("oris accelerate=TRUE weighted output sum = n (normalization preserved)", {
  fx <- load_stepstone()
  r  <- suppressWarnings(
    harvest(fx$df, fx$tgt, method = "oris", max_iterations = 50L, accelerate = TRUE)
  )
  w <- r$weights
  expect_equal(sum(w), nrow(fx$df), tolerance = 1e-6,
               label = "sum(weights) == n after SRAA")
})

## Test 8 — SOR coexists with SRAA: adapts on plain (non-AA) steps -----------

test_that("oris accelerate=TRUE allows SOR adaptation on plain SRAA steps", {
  fx <- load_stepstone()
  # Option B: SOR adaptation is re-enabled on plain (non-AA-accepted) SRAA
  # steps to dampen oscillating margins. AA-accepted steps still skip SOR
  # because their trajectory is non-monotone from extrapolation. As a result
  # sor_min_omega may drop below 1.0 when SRAA produces oscillating plain
  # steps. The contract: min_omega is in (0, 1] and finite.
  r  <- suppressWarnings(
    harvest(fx$df, fx$tgt, method = "oris", max_iterations = 100L, accelerate = TRUE)
  )
  min_omega <- attr(r, "result")$sor$min_omega
  expect_true(is.finite(min_omega),
              label = "sor_min_omega is finite under SRAA")
  expect_gt(min_omega, 0.0)
  expect_lte(min_omega, 1.0)
})

## Test 9 — Greedy scheduler + accelerate=TRUE downgrades to round_robin ----
# log goes via Rprintf (stdout), not R message conditions — use capture.output

test_that("oris greedy+accelerate downgrades to round_robin and logs at verbose=1", {
  set.seed(1L); n <- 200L
  df  <- data.frame(x = factor(sample(c("a","b"), n, TRUE)),
                    y = factor(sample(c("p","q"), n, TRUE)))
  tgt <- list(x = c(a=0.5, b=0.5), y = c(p=0.4, q=0.6))
  out <- capture.output(
    suppressWarnings(
      harvest(df, tgt, method = "oris", scheduler = "greedy",
              accelerate = TRUE, max_iterations = 20L, verbose = 1L)
    )
  )
  expect_true(any(grepl("round_robin", out)),
              label = "verbose=1 stdout contains 'round_robin'")
})
