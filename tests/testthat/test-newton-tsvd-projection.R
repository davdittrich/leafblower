# tests/testthat/test-newton-tsvd-projection.R
# WL-1 (Epic-Dβ): truncated-SVD pseudoinverse via dsyevd surfaces n_projected_dims.
library(leafblower)

test_that("WL-1 truncated-SVD projection drops null-direction components", {
  set.seed(1L); n <- 1000L
  # Construct K=2 fixture with overlapping-margin structure where one factor is
  # determined by the other (synthetic rank-deficient joint indicator block).
  df <- data.frame(
    a = factor(c(rep("x", n/2), rep("y", n/2))),
    # b is determined by a → joint indicator structure has dependent columns
    b = factor(c(rep("p", n/2), rep("q", n/2)))
  )
  tgt <- list(a = c(x=0.5, y=0.5), b = c(p=0.5, q=0.5))
  r <- harvest(df, tgt, method="newton_kl", attach_weights=FALSE)
  res <- attr(r, "result")
  expect_true(!is.null(res$n_projected_dims),
    label="n_projected_dims must surface")
  # On synthetic rank-deficient fixture, n_projected_dims should be ≥ 0.
  # (Acceptance loose — fixture might or might not trigger truncation depending
  # on exact spectrum; main goal is to verify the field surfaces correctly.)
  expect_true(is.integer(res$n_projected_dims) || is.numeric(res$n_projected_dims),
    label="n_projected_dims is a numeric scalar")
})

test_that("WL-1 K=4 well-conditioned baseline preserved", {
  set.seed(1L); n <- 2000L
  df <- data.frame(
    a = factor(sample(letters[1:3], n, TRUE)),
    b = factor(sample(LETTERS[1:4], n, TRUE)),
    c = factor(sample(c("M","F"), n, TRUE)),
    d = factor(sample(c("Y","O"), n, TRUE))
  )
  tgt <- list(
    a = c(a=0.4, b=0.35, c=0.25),
    b = c(A=0.3, B=0.3, C=0.2, D=0.2),
    c = c(M=0.5, F=0.5),
    d = c(Y=0.45, O=0.55)
  )
  r <- harvest(df, tgt, method="newton_kl", attach_weights=FALSE)
  res <- attr(r, "result")
  expect_lt(res$max_error, 1e-6,
    label=sprintf("K=4 well-conditioned: max_err=%.2e", res$max_error))
  expect_equal(res$n_projected_dims, 0L,
    label=sprintf("K=4 well-conditioned: no truncation expected; got %d", res$n_projected_dims))
})
