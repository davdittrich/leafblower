test_that("harvest(method='oris') returns algorithm='oris'", {
  set.seed(1L)
  n   <- 500L
  df  <- data.frame(x = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.4, b = 0.35, c = 0.25))
  r   <- harvest(df, tgt, method = "oris",
                 convergence = list(absolute = 1e-6),
                 attach_weights = FALSE)
  expect_equal(attr(r, "algorithm"), "oris")
})

test_that("harvest(method='oris_soft') returns algorithm='oris_soft'", {
  set.seed(2L)
  n   <- 500L
  df  <- data.frame(x = factor(sample(c("a", "b"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.6, b = 0.4))
  r   <- harvest(df, tgt, method = "oris_soft",
                 max_weight = 3, min_weight = 0,
                 convergence = list(absolute = 1e-6),
                 attach_weights = FALSE)
  expect_equal(attr(r, "algorithm"), "oris_soft")
})
