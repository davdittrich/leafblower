context("cell_table")

# Internal probe (test-only; returns M_cell, cell_of, n_per_cell)
probe <- function(group_ids_list, n) {
  .Call("C_leafblower_cell_table_probe", group_ids_list, as.integer(n),
        PACKAGE = "leafblower")
}

test_that("cell compression: all-identical observations produce 1 cell", {
  n <- 1000
  g1 <- rep(0L, n)
  g2 <- rep(0L, n)
  out <- probe(list(g1, g2), n)
  expect_equal(out$M_cell, 1L)
  expect_equal(out$n_per_cell, n)
  expect_true(all(out$cell_of == 0L))
})

test_that("cell compression: full cross-product populated", {
  # K=2 margins, 4 cats each, balanced assignment
  n <- 10000
  g1 <- rep(0:3, each = n/4)
  g2 <- rep(rep(0:3, each = n/16), 4)
  out <- probe(list(g1, g2), n)
  expect_equal(out$M_cell, 16L)  # 4 * 4
  expect_equal(sum(out$n_per_cell), n)
})

test_that("cell compression: sparse assignment gives M_cell < cross-product", {
  n <- 100
  g1 <- c(rep(0L, 50), rep(1L, 50))
  g2 <- c(rep(0L, 50), rep(1L, 50))  # perfectly correlated
  out <- probe(list(g1, g2), n)
  expect_equal(out$M_cell, 2L)  # only (0,0) and (1,1) populated
})

test_that("cell compression: NA entries (-1) treated as distinct category", {
  n <- 300
  g1 <- c(rep(0L, 100), rep(1L, 100), rep(-1L, 100))
  out <- probe(list(g1), n)
  expect_equal(out$M_cell, 3L)  # 0, 1, and NA buckets
})

test_that("cell compression: K > 64 rejected", {
  n <- 100
  many_margins <- replicate(65, rep(0L, n), simplify = FALSE)
  expect_error(probe(many_margins, n), regexp = "K.*64", ignore.case = TRUE)
})

test_that("cell compression: sort determinism (same input → same cell ids)", {
  n <- 500
  set.seed(1)
  g1 <- sample(0:3, n, replace = TRUE)
  g2 <- sample(0:2, n, replace = TRUE)
  out1 <- probe(list(g1, g2), n)
  out2 <- probe(list(g1, g2), n)
  expect_equal(out1$M_cell, out2$M_cell)
  expect_equal(out1$cell_of, out2$cell_of)
})

# eb79.4: TYPEOF guard on r_group_ids_list in the probe entry point. Pre-fix,
# a non-list group_ids goes straight into Rf_length()/VECTOR_ELT() on a
# non-VECSXP, which can segfault; isolate the check in a callr subprocess so
# a crash there doesn't abort the suite.
test_that("cell compression: non-list group_ids errors gracefully (no crash)", {
  n <- 100
  bad_group_ids <- as.integer(rep(0L, n))  # wrong TYPEOF: INTSXP vector, not a list
  run_wrong_typeof <- function(bad_group_ids, n) {
    library(leafblower)
    .Call("C_leafblower_cell_table_probe", bad_group_ids, as.integer(n),
          PACKAGE = "leafblower")
  }
  if (requireNamespace("callr", quietly = TRUE)) {
    expect_error(
      callr::r(run_wrong_typeof, args = list(bad_group_ids = bad_group_ids, n = n)),
      "group_ids must be a list"
    )
  } else {
    expect_error(run_wrong_typeof(bad_group_ids, n), "group_ids must be a list")
  }
})
