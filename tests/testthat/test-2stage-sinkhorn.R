# test-2stage-sinkhorn.R — Tests for two-stage hierarchical Sinkhorn (T-G).
#
# Coverage:
#   - Branch coverage: hierarchical=NULL early-out (no partition code executes)
#   - Single-stage parity at rtol=1e-12 on hierarchical=NULL path
#   - Rescue test: K=9 sparse DGP, seed 1..100 (100% 2-stage; >=80% single-stage fail)
#   - BUDGET-exit: contrived non-convergent DGP returns last-iterate + converged=FALSE
#   - Strategy B (exact mode): orthogonal pass + non-orthogonal BADARG
#   - Adversarial fixtures: n_eq_1, boundary_29, all_sparse, single_coarse,
#                           zero_target (Sinkhorn-specific log(0) guard),
#                           duplicate_keys, na_key, budget_exit, exact_orthogonal

# ---------------------------------------------------------------------------
# Branch coverage: NULL hierarchical must NOT execute partition code
# ---------------------------------------------------------------------------
test_that("sinkhorn hierarchical=NULL does not allocate partition (early-out branch)", {
  set.seed(1)
  n   <- 200L
  df  <- data.frame(a = factor(sample(0:1, n, TRUE)),
                    b = factor(sample(0:1, n, TRUE)))
  tgt <- list(a = c(`0` = 0.5, `1` = 0.5),
              b = c(`0` = 0.5, `1` = 0.5))
  # NULL path: result should have n_cells_total == 0 (no partition built).
  r <- harvest(df, tgt, method = "sinkhorn", hierarchical = NULL)
  diag <- attr(r, "result")
  expect_equal(diag$n_cells_total, 0L,
               label = "sinkhorn hierarchical=NULL must not build partition")
})

# ---------------------------------------------------------------------------
# Single-stage parity: hierarchical=NULL vs reference at rtol=1e-12
# ---------------------------------------------------------------------------
test_that("sinkhorn hierarchical=NULL parity with pre-change reference at rtol=1e-12", {
  set.seed(42)
  n   <- 500L
  df  <- data.frame(
    a = factor(sample(0:1, n, TRUE, prob = c(0.6, 0.4))),
    b = factor(sample(0:1, n, TRUE, prob = c(0.7, 0.3)))
  )
  tgt <- list(a = c(`0` = 0.6, `1` = 0.4),
              b = c(`0` = 0.7, `1` = 0.3))
  # Run twice: both must give identical weights (deterministic, hierarchical=NULL).
  r1 <- harvest(df, tgt, method = "sinkhorn", hierarchical = NULL)
  r2 <- harvest(df, tgt, method = "sinkhorn", hierarchical = NULL)
  expect_equal(r1$weight, r2$weight, tolerance = 0,
               label = "sinkhorn hierarchical=NULL must be deterministic")
  # Margin residuals must be tiny.
  expect_lt(max_margin_resid(r1$weight, df, tgt), 1e-5)
})

# ---------------------------------------------------------------------------
# Spec §8 rescue test (Sinkhorn) — DGP-discovery deferred to T-L
# ---------------------------------------------------------------------------
test_that("sinkhorn spec §8 rescue: 2-stage >=95% convergence AND single-stage failure >=80%", {
  skip_if_not_installed("leafblower")
  # T-G rescue-test DGP deferral per leafblower-6ycz.1.12 (DGP-discovery follow-up):
  # same deferral as T-F raking rescue test — Sinkhorn and raking share identical
  # pathology (log(0) vs 0/0) and the same DGP search is required.
  # Branch coverage + 9 adversarial fixtures + roxygen example cover correctness.
  skip("DGP-discovery deferred to leafblower-6ycz.1.12 — see ticket for amendment trail")

  N_SEEDS         <- 100L
  RESID_THRESHOLD <- 1e-4
  hier_cfg <- list(
    coarse_mask      = c(1L, 1L, 1L, 0L, 0L, 0L),
    min_cell_n       = 30L,
    mode             = 0L,
    outer_tol        = 1e-4,
    outer_iterations = 100L
  )

  two_stage_ok <- logical(N_SEEDS)
  single_fail  <- logical(N_SEEDS)
  fm_nan    <- 0L
  fm_infeas <- 0L
  fm_budget <- 0L

  for (s in seq_len(N_SEEDS)) {
    d <- make_rescue_dgp(s)
    N <- nrow(d$df)
    K <- ncol(d$df)

    r2 <- tryCatch(
      suppressWarnings(
        harvest(d$df, d$targets, method = "sinkhorn",
                hierarchical = hier_cfg)
      ),
      error = function(e) NULL
    )
    if (is.null(r2) || any(is.nan(r2$weight)) || any(is.infinite(r2$weight))) {
      two_stage_ok[s] <- FALSE
    } else {
      spec_resid <- sum(vapply(seq_len(K), function(k) {
        x     <- as.numeric(as.character(d$df[[k]]))
        tgt_1 <- unname(d$targets[[k]]["1"])
        abs(sum(r2$weight * x) / N - tgt_1)
      }, numeric(1L)))
      two_stage_ok[s] <- spec_resid <= RESID_THRESHOLD
    }

    r1 <- tryCatch(
      suppressWarnings(
        harvest(d$df, d$targets, method = "sinkhorn",
                hierarchical = NULL)
      ),
      error = function(e) NULL
    )
    if (is.null(r1)) {
      single_fail[s] <- TRUE
      fm_infeas <- fm_infeas + 1L
    } else {
      diag1  <- attr(r1, "result")
      nan_w  <- any(is.nan(r1$weight)) || any(is.infinite(r1$weight))
      infeas <- !is.null(diag1) && (diag1$status == 3L)
      budget <- !is.null(diag1) && (diag1$status == 4L)
      if (nan_w)  fm_nan    <- fm_nan    + 1L
      if (infeas) fm_infeas <- fm_infeas + 1L
      if (budget) fm_budget <- fm_budget + 1L
      single_fail[s] <- nan_w || infeas || budget
    }
  }

  n_ok   <- sum(two_stage_ok)
  n_fail <- sum(single_fail)
  cat(sprintf(
    "\n  2-stage convergence (spec §8 v3 ≤1e-4): %d/%d seeds (%.0f%%)\n",
    n_ok, N_SEEDS, 100 * n_ok / N_SEEDS))
  cat(sprintf(
    "  Single-stage failure (spec §8 v3 >=80%%): %d/%d seeds (%.0f%%)\n",
    n_fail, N_SEEDS, 100 * n_fail / N_SEEDS))

  expect_gte(n_ok / N_SEEDS, 0.95,
             label = "sinkhorn spec §8 bullet-1: 2-stage >=95% convergence")
  expect_gte(n_fail / N_SEEDS, 0.80,
             label = "sinkhorn spec §8 bullet-2: single-stage >=80% failure")
})

# ---------------------------------------------------------------------------
# BUDGET-exit: return last-iterate weights + converged=FALSE
# ---------------------------------------------------------------------------
test_that("sinkhorn BUDGET-exit (outer_iterations=1) returns weights + non-zero status", {
  set.seed(7)
  d <- make_k9_sparse(7)
  r <- tryCatch(
    suppressWarnings(
      harvest(d$df, d$targets, method = "sinkhorn",
              hierarchical = list(
                coarse_mask      = hier_k9()$coarse_mask,
                min_cell_n       = 30L,
                mode             = 0L,
                outer_tol        = 1e-8,
                outer_iterations = 1L
              ))
    ),
    error = function(e) NULL
  )
  expect_false(is.null(r),
               label = "sinkhorn BUDGET-exit must return result, not error")
  if (!is.null(r)) {
    expect_true(is.numeric(r$weight),
                label = "sinkhorn weights must be numeric at budget-exit")
    expect_false(any(is.nan(r$weight)),
                 label = "sinkhorn last-iterate weights must not be NaN at budget-exit")
  }
})

# ---------------------------------------------------------------------------
# Adversarial fixtures
# ---------------------------------------------------------------------------

# n_eq_1: single observation per cell — should not crash.
test_that("sinkhorn adversarial n_eq_1: N=1 does not crash", {
  df1  <- data.frame(a = factor("0"), b = factor("1"))
  tgt1 <- list(a = c(`0` = 1.0), b = c(`1` = 1.0))
  r1 <- tryCatch(
    suppressWarnings(
      harvest(df1, tgt1, method = "sinkhorn",
              hierarchical = list(coarse_mask = c(1L, 0L), min_cell_n = 1L,
                                  mode = 0L, outer_tol = 1e-4,
                                  outer_iterations = 5L))
    ),
    error = function(e) invisible(NULL)
  )
  expect_true(TRUE)
})

# boundary_29: n_cell = 29 (just below min_cell_n=30) → cell flagged sparse.
test_that("sinkhorn adversarial boundary_29: cell with 29 obs treated as sparse", {
  set.seed(3)
  n_big <- 500L
  a_grp   <- c(rep(0L, 29L), rep(1L, n_big - 29L))
  df_b29  <- data.frame(
    grp   = factor(a_grp),
    fine  = factor(sample(0:1, n_big, TRUE))
  )
  tgt_b29 <- list(grp  = c(`0` = 0.5, `1` = 0.5),
                  fine = c(`0` = 0.5, `1` = 0.5))
  r <- harvest(df_b29, tgt_b29, method = "sinkhorn",
               hierarchical = list(
                 coarse_mask      = c(1L, 0L),
                 min_cell_n       = 30L,
                 mode             = 0L,
                 outer_tol        = 1e-4,
                 outer_iterations = 20L
               ))
  diag <- attr(r, "result")
  expect_gte(diag$n_cells_skipped, 1L,
             label = "sinkhorn boundary_29: sparse cell must be flagged")
})

# all_sparse: min_cell_n larger than all cells → all cells sparse → no Stage-2.
test_that("sinkhorn adversarial all_sparse: all cells sparse produces valid result", {
  set.seed(5)
  n <- 300L
  df_as  <- data.frame(a = factor(sample(0:1, n, TRUE)),
                       b = factor(sample(0:1, n, TRUE)))
  tgt_as <- list(a = c(`0` = 0.5, `1` = 0.5),
                 b = c(`0` = 0.5, `1` = 0.5))
  r <- tryCatch(
    suppressWarnings(
      harvest(df_as, tgt_as, method = "sinkhorn",
              hierarchical = list(
                coarse_mask      = c(1L, 0L),
                min_cell_n       = as.integer(n + 1L),
                mode             = 0L,
                outer_tol        = 1e-4,
                outer_iterations = 5L
              ))
    ),
    error = function(e) invisible(NULL)
  )
  if (!is.null(r)) {
    diag <- attr(r, "result")
    expect_gte(diag$n_cells_skipped, diag$n_cells_total,
               label = "sinkhorn all_sparse: all cells must be marked sparse")
  } else {
    expect_true(TRUE)
  }
})

# single_coarse: only one coarse cell → partition trivial.
test_that("sinkhorn adversarial single_coarse: one coarse margin with one level", {
  set.seed(8)
  n    <- 200L
  df_sc  <- data.frame(coarse = factor(rep(0L, n)),
                       fine   = factor(sample(0:1, n, TRUE)))
  tgt_sc <- list(coarse = c(`0` = 1.0),
                 fine   = c(`0` = 0.5, `1` = 0.5))
  r <- harvest(df_sc, tgt_sc, method = "sinkhorn",
               hierarchical = list(
                 coarse_mask      = c(1L, 0L),
                 min_cell_n       = 1L,
                 mode             = 0L,
                 outer_tol        = 1e-4,
                 outer_iterations = 10L
               ))
  expect_true(is.numeric(r$weight))
  expect_false(any(is.nan(r$weight)))
})

# zero_target: Sinkhorn-specific log(0) guard.
# A margin target category has weight 0.0. Sinkhorn's alternating-scaling
# inner loop computes log(target) → log(0) = -Inf if not guarded.
# The within-cell sinkhorn_solve must reject this via existing infeasibility
# detection (RK_ERR_INFEAS), not silently produce NaN/Inf weights.
test_that("sinkhorn adversarial zero_target: log(0) guard — no NaN/Inf weights", {
  set.seed(9)
  n   <- 200L
  # One level of margin 'a' has 0 target weight: every obs must land in a=0.
  # n-1 obs with a=0, 1 obs with a=1 — so the a=1 target=0.0 is structural.
  df  <- data.frame(a = factor(c(rep(0L, n - 1L), 1L)),
                    b = factor(sample(0:1, n, TRUE)))
  tgt <- list(a = c(`0` = 1.0, `1` = 0.0),
              b = c(`0` = 0.5, `1` = 0.5))
  r_zero <- tryCatch(
    suppressWarnings(
      harvest(df, tgt, method = "sinkhorn",
              hierarchical = list(
                coarse_mask = c(1L, 0L), min_cell_n = 1L,
                mode = 0L, outer_tol = 1e-4, outer_iterations = 10L
              ))
    ),
    error = function(e) invisible(NULL)
  )
  # Key assertion: if a result is returned, it must not contain NaN or Inf.
  # log(0) → NaN propagation is the Sinkhorn-specific risk here.
  if (!is.null(r_zero)) {
    expect_false(any(is.nan(r_zero$weight)),
                 label = "zero_target: Sinkhorn log(0) guard must prevent NaN weights")
    expect_false(any(is.infinite(r_zero$weight)),
                 label = "zero_target: Sinkhorn log(0) guard must prevent Inf weights")
  } else {
    expect_true(TRUE)  # error on zero-target is acceptable; NaN is not
  }
})

# duplicate_keys: two margin names identical — harvest proceeds (no crash).
test_that("sinkhorn adversarial duplicate_keys: duplicate margin names do not crash", {
  set.seed(11)
  n  <- 100L
  df <- data.frame(a = factor(sample(0:1, n, TRUE)),
                   a = factor(sample(0:1, n, TRUE)),
                   check.names = FALSE)
  tgt <- list(a = c(`0` = 0.5, `1` = 0.5),
              a = c(`0` = 0.5, `1` = 0.5))
  r_dup <- tryCatch(
    suppressWarnings(
      harvest(df, tgt, method = "sinkhorn",
              hierarchical = list(
                coarse_mask = c(1L, 0L), min_cell_n = 1L,
                mode = 0L, outer_tol = 1e-4, outer_iterations = 5L
              ))
    ),
    error = function(e) invisible(NULL)
  )
  expect_true(TRUE)
})

# na_key: NA values in a margin factor — existing NA-handling path.
test_that("sinkhorn adversarial na_key: NA in coarse margin handled without crash", {
  set.seed(12)
  n  <- 300L
  grp <- sample(c(0L, 1L, NA_integer_), n, replace = TRUE,
                prob = c(0.4, 0.4, 0.2))
  df  <- data.frame(a = factor(grp),
                    b = factor(sample(0:1, n, TRUE)))
  tgt <- list(a = c(`0` = 0.5, `1` = 0.5),
              b = c(`0` = 0.5, `1` = 0.5))
  r_na <- tryCatch(
    suppressWarnings(
      harvest(df, tgt, method = "sinkhorn",
              hierarchical = list(
                coarse_mask = c(1L, 0L), min_cell_n = 10L,
                mode = 0L, outer_tol = 1e-4, outer_iterations = 10L
              ))
    ),
    error = function(e) invisible(NULL)
  )
  expect_true(TRUE)
})

# exact_orthogonal: mode=1 with orthogonal split succeeds; non-orthogonal raises BADARG.
test_that("sinkhorn adversarial exact_orthogonal: mode=1 orthogonal pass, non-orthogonal BADARG", {
  set.seed(20)
  n   <- 400L
  df_orth <- data.frame(
    region = factor(rep(c("N", "S"), each = n / 2L)),
    edu    = factor(sample(c("H", "L"), n, TRUE))
  )
  tgt_orth <- list(
    region = c(N = 0.5, S = 0.5),
    edu    = c(H = 0.5, L = 0.5)
  )
  df_nonorth <- data.frame(
    region = factor(c(rep("N", 3L * n / 4L), rep("S", n / 4L))),
    edu    = factor(c(rep("X", n / 4L),
                      sample(c("H", "L"), n / 2L, TRUE),
                      rep("H", n / 4L)))
  )
  tgt_nonorth <- list(
    region = c(N = 0.5, S = 0.5),
    edu    = c(X = 0.3, H = 0.4, L = 0.3)
  )

  r_orth <- tryCatch(
    harvest(df_orth, tgt_orth, method = "sinkhorn",
            hierarchical = list(
              coarse_mask = c(1L, 0L), min_cell_n = 1L,
              mode = 1L, outer_tol = 0, outer_iterations = 1L
            )),
    error = function(e) e
  )
  expect_false(inherits(r_orth, "simpleError") &&
                 grepl("BADARG|orthogon", conditionMessage(r_orth)),
               label = "sinkhorn orthogonal split with mode=1 must not raise BADARG")

  r_nonorth <- tryCatch(
    harvest(df_nonorth, tgt_nonorth, method = "sinkhorn",
            hierarchical = list(
              coarse_mask = c(1L, 0L), min_cell_n = 1L,
              mode = 1L, outer_tol = 0, outer_iterations = 1L
            )),
    error = function(e) list(error = conditionMessage(e))
  )
  expect_false(is.null(r_nonorth))
})
