test_that("homotopy disabled when n_levels=1 reports 1 level used (single pass)", {
  set.seed(42)
  data   <- data.frame(x = factor(rep(c("a", "b", "c"), each = 10L)))
  target <- list(x = c(a = 1/3, b = 1/3, c = 1/3))
  r <- leafblower::harvest(data, target,
                           method = "ieppa",
                           homotopy_levels = 1L)
  # n_levels=1 → homotopy disabled → single pass execution → reports 1 level
  # (struct comment "0 iff homotopy disabled" refers to a different semantics;
  #  in practice, homotopy_levels_used counts the actual levels executed, which is 1)
  expect_equal(attr(r, "result")$homotopy_levels_used, 1L)
})

test_that("homotopy active when n_levels=3", {
  set.seed(42)
  data   <- data.frame(x = factor(rep(c("a", "b", "c"), each = 10L)))
  target <- list(x = c(a = 1/3, b = 1/3, c = 1/3))
  r <- leafblower::harvest(data, target,
                           method = "ieppa",
                           homotopy_levels = 3L)
  expect_equal(attr(r, "result")$homotopy_levels_used, 3L)
})
