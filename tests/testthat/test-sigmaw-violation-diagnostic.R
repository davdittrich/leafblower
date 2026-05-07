## T-P: enforce_sigmaw_eq_n_diag probe tests
## Verifies the new diagnostic overload (leafblower-6ycz.1.16).

diag_probe <- function(weights, N) {
    .Call("C_hier_sigmaw_diag_probe", as.double(weights), as.integer(N))
}

test_that("diag probe: weights sum exactly to N → passed=TRUE, dev≈0", {
    w <- c(1.0, 1.0, 1.0)
    res <- diag_probe(w, 3L)
    expect_true(res$passed)
    expect_equal(res$dev, 0.0, tolerance = 1e-15)
})

test_that("diag probe: weights do not sum to N → passed=FALSE, correct dev", {
    # sum = 3.5, N = 3, dev = 0.5
    w <- c(1.0, 1.0, 1.5)
    res <- diag_probe(w, 3L)
    expect_false(res$passed)
    expect_equal(res$dev, 0.5, tolerance = 1e-15)
})

test_that("diag probe: small violation below tolerance → passed=TRUE", {
    # tol = 3 * 1e-12 = 3e-12; dev = 1e-13 < tol
    N <- 3L
    w <- c(1.0, 1.0, 1.0 + 1e-13)
    res <- diag_probe(w, N)
    expect_true(res$passed)
    expect_lt(res$dev, N * 1e-12)
})

test_that("diag probe: violation just above tolerance → passed=FALSE", {
    # tol = 3 * 1e-12; dev = 1e-11 > tol
    N <- 3L
    w <- c(1.0, 1.0, 1.0 + 1e-11)
    res <- diag_probe(w, N)
    expect_false(res$passed)
    expect_gt(res$dev, N * 1e-12)
})

test_that("diag probe: N <= 0 → passed=FALSE, dev=Inf", {
    w <- c(1.0, 1.0, 1.0)
    res <- diag_probe(w, 0L)
    expect_false(res$passed)
    expect_true(is.infinite(res$dev))
})

test_that("diag probe: result list has names 'passed' and 'dev'", {
    w <- c(1.0, 2.0)
    res <- diag_probe(w, 2L)
    expect_named(res, c("passed", "dev"))
})
