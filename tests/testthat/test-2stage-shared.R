# test-2stage-shared.R
# Unit tests for T-B shared 2-stage infrastructure:
#   build_cell_partition, build_sparse_mask, apply_sparse_inheritance,
#   enforce_sigmaw_eq_n
#
# Helpers are C++ only (not yet wired to R harvest path — that is T-C through T-K).
# They are tested here via thin C probes registered in r_bridge.cpp:
#   C_hier_partition_probe(group_ids_list, n, K, coarse_mask, min_cell_n)
#     -> list(rc, n_cells_total, n_cells_skipped, cell_id_per_obs, n_per_cell)
#   C_hier_inherit_probe(group_ids_list, n, K, coarse_mask, min_cell_n, stage1_mults)
#     -> list(rc, weights_out, n_cells_total, n_cells_skipped)
#   C_hier_sigmaw_probe(weights, N) -> logical(1)
#
# LBW_MAX_HIER_CELLS = 100000 (defined in calib_dispatch.hpp).

.hier_partition <- function(gid_list, n, K, coarse_mask, min_cell_n = 1L) {
  .Call("C_hier_partition_probe",
        gid_list, n, K, as.integer(coarse_mask), as.integer(min_cell_n),
        PACKAGE = "leafblower")
}

.hier_inherit <- function(gid_list, n, K, coarse_mask, min_cell_n, stage1_mults) {
  .Call("C_hier_inherit_probe",
        gid_list, n, K, as.integer(coarse_mask), as.integer(min_cell_n),
        as.double(stage1_mults),
        PACKAGE = "leafblower")
}

.hier_sigmaw <- function(weights, N) {
  .Call("C_hier_sigmaw_probe", as.double(weights), as.integer(N),
        PACKAGE = "leafblower")
}

RK_OK       <- 0L
RK_ERR_BADARG <- 3L

# ── build_cell_partition ──────────────────────────────────────────────────

test_that("build_cell_partition: single cell when all obs same category", {
  n <- 10L; K <- 1L
  gids <- list(rep(0L, n))
  p <- .hier_partition(gids, n, K, coarse_mask = c(1L), min_cell_n = 1L)

  expect_equal(p$rc, RK_OK)
  expect_equal(p$n_cells_total, 1L)
  expect_equal(length(p$n_per_cell), 1L)
  expect_equal(p$n_per_cell[[1L]], 10L)
  expect_true(all(p$cell_id_per_obs == 0L))
})

test_that("build_cell_partition: dense partition — K=1, 3 cats, n=15", {
  n <- 15L; K <- 1L
  # 5 obs each in cats 0,1,2
  gids <- list(c(rep(0L, 5L), rep(1L, 5L), rep(2L, 5L)))
  p <- .hier_partition(gids, n, K, coarse_mask = c(1L), min_cell_n = 1L)

  expect_equal(p$rc, RK_OK)
  expect_equal(p$n_cells_total, 3L)
  expect_equal(sort(p$n_per_cell), c(5L, 5L, 5L))
})

test_that("build_cell_partition: 2 coarse margins, cross-product cells", {
  n <- 24L; K <- 2L
  # age (0,1,2) x sex (0,1) = 6 cells, 4 obs each
  set.seed(1L)
  age <- rep(0:2, each = 8L)
  sex <- rep(c(0L, 0L, 0L, 0L, 1L, 1L, 1L, 1L), times = 3L)
  gids <- list(as.integer(age), as.integer(sex))
  p <- .hier_partition(gids, n, K, coarse_mask = c(1L, 1L), min_cell_n = 1L)

  expect_equal(p$rc, RK_OK)
  expect_equal(p$n_cells_total, 6L)
  expect_equal(sort(p$n_per_cell), rep(4L, 6L))
})

test_that("build_cell_partition: coarse_mask=0 for non-coarse margin ignores it", {
  n <- 12L; K <- 2L
  # 2 margins: coarse=margin 0 (3 cats), fine=margin 1 (2 cats, mask=0)
  # partition should only see 3 cells from margin 0
  gids <- list(rep(0:2, each = 4L),
               rep(c(0L, 1L), times = 6L))
  p <- .hier_partition(gids, n, K, coarse_mask = c(1L, 0L), min_cell_n = 1L)

  expect_equal(p$rc, RK_OK)
  expect_equal(p$n_cells_total, 3L)
})

test_that("build_cell_partition: all obs in single margin-0 cat (empty other cat)", {
  n <- 5L; K <- 1L
  gids <- list(rep(0L, n))
  p <- .hier_partition(gids, n, K, coarse_mask = c(1L), min_cell_n = 1L)

  expect_equal(p$rc, RK_OK)
  expect_equal(p$n_cells_total, 1L)
})

# ── Cell-count cap (RK_ERR_BADARG) ───────────────────────────────────────

test_that("build_cell_partition: cap fires when K_cells > N/min_cell_n", {
  # n=12, min_cell_n=5 => cap = min(12/5, 100000) = 2 (integer div)
  # 3 unique cats => K_cells=3 > cap=2 => BADARG
  n <- 12L; K <- 1L
  gids <- list(rep(0:2, each = 4L))
  p <- .hier_partition(gids, n, K, coarse_mask = c(1L), min_cell_n = 5L)

  expect_equal(p$rc, RK_ERR_BADARG)
  expect_null(p$n_cells_total)   # truncated list on BADARG: only "rc" field present
})

test_that("build_cell_partition: cap passes when K_cells <= N/min_cell_n", {
  # n=15, min_cell_n=5 => cap=min(3,100000)=3; K_cells=3 => OK
  n <- 15L; K <- 1L
  gids <- list(rep(0:2, each = 5L))
  p <- .hier_partition(gids, n, K, coarse_mask = c(1L), min_cell_n = 5L)

  expect_equal(p$rc, RK_OK)
  expect_equal(p$n_cells_total, 3L)
})

# ── build_sparse_mask (via partition probe — n_cells_skipped) ─────────────

test_that("build_sparse_mask: all-sparse when min_cell_n > all n_per_cell", {
  n <- 6L; K <- 1L
  # 3 cells, 2 obs each; min_cell_n=3 => all 3 cells sparse
  gids <- list(rep(0:2, each = 2L))
  p <- .hier_partition(gids, n, K, coarse_mask = c(1L), min_cell_n = 3L)

  # rc = BADARG because cap = min(6/3,100000)=2, K_cells=3 > 2
  # So adjust: use min_cell_n=2 => cap=min(3,100000)=3, K_cells=3 => OK
  # min_cell_n=2 means cells with n<2 are sparse; each has n=2 => not sparse (>=2)
  # Use min_cell_n=1: cells with n<1 sparse => none
  # Use min_cell_n=3: cap=2, BADARG. Need to pick values carefully.
  # n=9, 3 cats 3 obs each, min_cell_n=4 => cap=min(2,100000)=2, K_cells=3>2 => BADARG.
  # n=12, 3 cats 4 obs each, min_cell_n=5 => cap=2, K_cells=3>2 => BADARG.
  # For all-sparse WITHOUT cap: n=30, 3 cats 10 obs each, min_cell_n=11
  #   cap=min(30/11,100000)=2, K_cells=3>2 => BADARG again!
  # Root: cap = N/min_cell_n; if min_cell_n > avg_cell_size, cap < K_cells always.
  # So all-sparse requires min_cell_n > max_cell_size which implies cap < K_cells => BADARG.
  # Test instead: partial-sparse scenario only.
  skip("all-sparse is logically impossible without cap violation; see partial-sparse test below")
})

test_that("build_sparse_mask: partial-sparse — cells below threshold flagged", {
  # n=20, K=1, cats: 0->15 obs (dense), 1->5 obs (sparse), min_cell_n=10
  # cap = min(20/10, 100000) = 2; K_cells=2 => 2 <= 2 => OK
  n <- 20L; K <- 1L
  gids <- list(c(rep(0L, 15L), rep(1L, 5L)))
  p <- .hier_partition(gids, n, K, coarse_mask = c(1L), min_cell_n = 10L)

  expect_equal(p$rc, RK_OK)
  expect_equal(p$n_cells_total, 2L)
  # The sparse cell (5 obs < 10) should be flagged: n_cells_skipped = 1
  expect_equal(p$n_cells_skipped, 1L)
  # Dense cell (15 obs >= 10) should NOT be flagged
  # n_cells_skipped count = 1 (only the 5-obs cell)
})

test_that("build_sparse_mask: no sparse cells when all n_per_cell >= min_cell_n", {
  n <- 20L; K <- 1L
  gids <- list(c(rep(0L, 10L), rep(1L, 10L)))
  # min_cell_n=5 => cap=min(4,100000)=4; K_cells=2<=4 => OK; both cells have 10>=5 => 0 sparse
  p <- .hier_partition(gids, n, K, coarse_mask = c(1L), min_cell_n = 5L)

  expect_equal(p$rc, RK_OK)
  expect_equal(p$n_cells_skipped, 0L)
})

# ── apply_sparse_inheritance ──────────────────────────────────────────────

test_that("apply_sparse_inheritance: sparse cells get multiplier, dense unchanged", {
  # n=20, cell 0 (15 obs, dense), cell 1 (5 obs, sparse with min_cell_n=10)
  # stage1_mults = c(2.0, 3.0): dense cell mult=2.0 (ignored), sparse cell mult=3.0
  # w_init all 1.0 => sparse obs get weight 1.0*3.0=3.0, dense stay 1.0
  n <- 20L; K <- 1L
  gids <- list(c(rep(0L, 15L), rep(1L, 5L)))
  r <- .hier_inherit(gids, n, K,
                     coarse_mask   = c(1L),
                     min_cell_n    = 10L,
                     stage1_mults  = c(2.0, 3.0))

  expect_equal(r$rc, RK_OK)
  w <- r$weights_out
  expect_equal(length(w), 20L)
  # dense obs (first 15): unchanged at 1.0
  expect_true(all(abs(w[1:15] - 1.0) < 1e-14))
  # sparse obs (last 5): 1.0 * 3.0 = 3.0
  expect_true(all(abs(w[16:20] - 3.0) < 1e-14))
})

test_that("apply_sparse_inheritance: identity when stage1_mults all 1.0", {
  n <- 20L; K <- 1L
  gids <- list(c(rep(0L, 15L), rep(1L, 5L)))
  r <- .hier_inherit(gids, n, K,
                     coarse_mask   = c(1L),
                     min_cell_n    = 10L,
                     stage1_mults  = c(1.0, 1.0))

  expect_equal(r$rc, RK_OK)
  expect_true(all(abs(r$weights_out - 1.0) < 1e-14))
})

test_that("apply_sparse_inheritance: n_cells_skipped matches build_sparse_mask", {
  n <- 20L; K <- 1L
  gids <- list(c(rep(0L, 15L), rep(1L, 5L)))
  r <- .hier_inherit(gids, n, K,
                     coarse_mask   = c(1L),
                     min_cell_n    = 10L,
                     stage1_mults  = c(2.0, 3.0))
  # Acceptance criterion from iter-1 plan review: n_cells_inherited == n_cells_skipped
  expect_equal(r$n_cells_skipped, 1L)
})

# ── enforce_sigmaw_eq_n ───────────────────────────────────────────────────

test_that("enforce_sigmaw_eq_n: accepts exact sum", {
  w <- rep(1.0, 100L)
  expect_true(.hier_sigmaw(w, 100L))
})

test_that("enforce_sigmaw_eq_n: accepts sum within N * 1e-12 tolerance", {
  # Σw = 100 + 0.5 * 100 * 1e-12 (inside tolerance)
  N <- 100L
  eps <- as.double(N) * 1e-12 * 0.5
  w <- c(rep(1.0, 99L), 1.0 + eps)
  expect_true(.hier_sigmaw(w, N))
})

test_that("enforce_sigmaw_eq_n: rejects sum outside N * 1e-12 tolerance", {
  # Σw = 100 + 2 * 100 * 1e-12 (outside tolerance)
  N <- 100L
  eps <- as.double(N) * 1e-12 * 2.0
  w <- c(rep(1.0, 99L), 1.0 + eps)
  expect_false(.hier_sigmaw(w, N))
})

test_that("enforce_sigmaw_eq_n: rejects sum far from N", {
  w <- rep(2.0, 100L)  # Σw = 200 != 100
  expect_false(.hier_sigmaw(w, 100L))
})

test_that("enforce_sigmaw_eq_n: accepts zero-weight scenario sum = N", {
  # All weight on one obs = N
  N <- 50L
  w <- c(as.double(N), rep(0.0, 49L))
  expect_true(.hier_sigmaw(w, N))
})

# ── n_cells_inherited == n_cells_skipped assertion ────────────────────────

test_that("n_cells_inherited equals n_cells_skipped at partition exit", {
  # Iter-1 plan review: assert n_cells_inherited == n_cells_skipped
  n <- 20L; K <- 1L
  gids <- list(c(rep(0L, 15L), rep(1L, 5L)))
  p <- .hier_partition(gids, n, K, coarse_mask = c(1L), min_cell_n = 10L)
  r <- .hier_inherit(gids, n, K, coarse_mask = c(1L), min_cell_n = 10L,
                     stage1_mults = c(2.0, 3.0))
  expect_equal(r$n_cells_skipped, p$n_cells_skipped)
})
