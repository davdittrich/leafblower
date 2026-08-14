# These tests verify that the ALMConfig substruct refactor (st.alm.mu,
# st.alm.capacity_mu, st.alm.lambda) did not break solver behaviour.
# No behavioral change is expected — only access-path renamed in C++.

test_that("oris_soft converges after ALMConfig grouping", {
  set.seed(7)
  n  <- 200L
  df <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  r <- harvest(df, tgt, method = "oris_soft", max_iterations = 200L)
  expect_lte(attr(r, "result")$max_error, 1e-3)
})

test_that("oris converges (alm inactive) after ALMConfig grouping", {
  set.seed(7)
  n  <- 200L
  df <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  r <- harvest(df, tgt, method = "oris", max_iterations = 200L)
  expect_lte(attr(r, "result")$max_error, 1e-3)
})
