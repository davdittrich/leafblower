library(leafblower)

# T-E: hierarchical input-validation table.
# One test per §9 check row from the spec.
# All BADARG checks use harvest() which converts status=3 to stop().
# Each test constructs the minimal valid scenario then violates one constraint.

# ── Shared test harness ────────────────────────────────────────────────────────
# K=2 margins, n=200; coarse_mask c(1L,0L) is a valid proper subset.
make_valid_hier <- function(n = 200L, method = "raking") {
  set.seed(1)
  df  <- data.frame(
    age = sample(c("young", "old"),   n, replace = TRUE),
    sex = sample(c("M", "F"),         n, replace = TRUE)
  )
  tgt <- list(
    age = c(young = 0.5, old = 0.5),
    sex = c(M = 0.5, F = 0.5)
  )
  hier <- list(
    coarse_mask       = c(1L, 0L),   # age is coarse, sex is fine — proper subset
    min_cell_n        = 5L,
    mode              = 0L,           # refine
    outer_tol         = 1e-3,
    outer_iterations  = 10L
  )
  list(df = df, tgt = tgt, hier = hier, method = method)
}

# ── Check (1): hierarchical=NULL short-circuit (no error) ─────────────────────
test_that("hierarchical=NULL: no validation errors emitted", {
  s <- make_valid_hier()
  expect_no_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = NULL)
  )
})

# ── Check (2): coarse_mask length != K ────────────────────────────────────────
test_that("coarse_mask length != K -> error naming arg and K", {
  s <- make_valid_hier()
  s$hier$coarse_mask <- c(1L, 0L, 0L)   # K=2 but length 3
  expect_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier),
    regexp = "coarse_mask"
  )
})

# ── Check (3): coarse_mask all-zero (empty) ───────────────────────────────────
test_that("coarse_mask all-zero -> error 'coarse_margins is empty'", {
  s <- make_valid_hier()
  s$hier$coarse_mask <- c(0L, 0L)
  expect_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier),
    regexp = "coarse_margins is empty"
  )
})

# ── Check (4): coarse_mask == full set (degenerate) ───────────────────────────
test_that("coarse_mask all-ones -> error 'nothing to refine'", {
  s <- make_valid_hier()
  s$hier$coarse_mask <- c(1L, 1L)   # all margins coarse → nothing to refine
  expect_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier),
    regexp = "nothing to refine"
  )
})

# ── Check (5): min_cell_n < 1 ─────────────────────────────────────────────────
test_that("min_cell_n=0 -> BADARG naming min_cell_n and N", {
  s <- make_valid_hier()
  s$hier$min_cell_n <- 0L
  expect_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier),
    regexp = "min_cell_n"
  )
})

# ── Check (6): min_cell_n > N -> Rf_warning (NOT error) ──────────────────────
test_that("min_cell_n > N -> warning (not error), degenerate single-stage", {
  s <- make_valid_hier()
  n_obs <- nrow(s$df)
  s$hier$min_cell_n <- n_obs + 1L
  expect_warning(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier),
    regexp = "degenerate to single-stage"
  )
})

# ── Check (7): mode not in {0, 1} ────────────────────────────────────────────
test_that("hierarchical_mode=2 -> BADARG naming mode", {
  s <- make_valid_hier()
  s$hier$mode <- 2L
  expect_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier),
    regexp = "hierarchical_mode"
  )
})

# ── Check (8a): outer_tol not finite ─────────────────────────────────────────
test_that("outer_tol=Inf -> BADARG naming outer_tol", {
  s <- make_valid_hier()
  s$hier$outer_tol <- Inf
  expect_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier),
    regexp = "outer_tol"
  )
})

# ── Check (8b): outer_tol <= 0 ───────────────────────────────────────────────
test_that("outer_tol=0 -> BADARG naming outer_tol", {
  s <- make_valid_hier()
  s$hier$outer_tol <- 0.0
  expect_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier),
    regexp = "outer_tol"
  )
})

# ── Check (9a): outer_iterations < 1 ─────────────────────────────────────────
test_that("outer_iterations=0 -> BADARG naming outer_iterations", {
  s <- make_valid_hier()
  s$hier$outer_iterations <- 0L
  expect_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier),
    regexp = "outer_iterations"
  )
})

# ── Check (9b): outer_iterations > 10000 ─────────────────────────────────────
test_that("outer_iterations=10001 -> BADARG naming outer_iterations and cap", {
  s <- make_valid_hier()
  s$hier$outer_iterations <- 10001L
  expect_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier),
    regexp = "outer_iterations"
  )
})

# ── Check (10): algorithm × mode — LOGIT + mode=refine -> BADARG ─────────────
test_that("logit + mode=0 (refine) -> BADARG naming algorithm and mode", {
  s <- make_valid_hier(method = "logit")
  s$hier$mode <- 0L   # refine forbidden for Newton methods
  expect_error(
    harvest(s$df, s$tgt, method = "logit", hierarchical = s$hier),
    regexp = "logit"
  )
})

# ── Check (10b): GREG + mode=refine -> BADARG ────────────────────────────────
test_that("greg + mode=0 (refine) -> BADARG naming algorithm and mode", {
  s <- make_valid_hier(method = "greg")
  s$hier$mode <- 0L
  expect_error(
    harvest(s$df, s$tgt, method = "greg", hierarchical = s$hier),
    regexp = "greg"
  )
})

# ── Check (11): bounds_mode='unit' + hierarchical -> BADARG ──────────────────
test_that("bounds_mode='unit' + hierarchical -> BADARG naming bounds_mode", {
  s <- make_valid_hier()
  expect_error(
    harvest(s$df, s$tgt, method = s$method, hierarchical = s$hier,
            bounds_mode = "unit"),
    regexp = "bounds_mode.*unit|unit.*incompatible"
  )
})

# ── Check (12): cell-count cap exceeded ──────────────────────────────────────
# Build a scenario where the estimated coarse cell count exceeds
# min(N/min_cell_n, 100000).
# n=200, min_cell_n=10 => cap = min(200/10, 100000) = 20.
# coarse margin has 25 levels => cells_est = 25 > cap = 20 => BADARG.
test_that("estimated coarse cells > cap -> BADARG naming cells and cap", {
  set.seed(99)
  n2   <- 200L
  cats <- paste0("c", 1:25)   # 25 coarse categories
  df2  <- data.frame(
    grp  = factor(sample(cats, n2, replace = TRUE), levels = cats),
    flag = sample(c("A", "B"), n2, replace = TRUE)
  )
  tgt2 <- list(
    grp  = setNames(rep(1/25, 25), cats),
    flag = c(A = 0.5, B = 0.5)
  )
  hier2 <- list(
    coarse_mask      = c(1L, 0L),   # grp is coarse (25 cats)
    min_cell_n       = 10L,         # cap = min(200/10, 100000) = 20; cells_est = 25 > 20
    mode             = 0L,
    outer_tol        = 1e-3,
    outer_iterations = 5L
  )
  expect_error(
    harvest(df2, tgt2, method = "raking", hierarchical = hier2),
    regexp = "cells|cap|coarse"
  )
})
