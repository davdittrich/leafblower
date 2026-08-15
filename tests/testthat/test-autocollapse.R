# Under testthat 3e, expect_warning() only muffles the ONE warning matching
# its regexp; any other warning from the same call (here harvest()'s
# incidental sparse-category diagnostic) propagates -- edition 2 used to
# swallow all of them (D-11, phase 01-04). Muffle ONLY that named incidental
# pattern rather than every warning, so edition 3's stricter unmatched-warning
# signal still catches an unrelated regression in these calls.
.muffle_sparse_category_diag <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("sparse categories detected", conditionMessage(w), fixed = TRUE))
      invokeRestart("muffleWarning")
  })
}

test_that("auto_collapse handles rare categories without crashing (eb79.17)", {
  df <- data.frame(x = factor(c(rep("a", 100), rep("b", 95), rep("z", 5))))
  res <- harvest(df, list(x = c(a = .45, b = .45, z = .10)),
                 auto_collapse = TRUE, max_weight = 5)

  expect_equal(attr(res, "result")$status, 0L)   # 0 = RK_OK
  expect_true(all(is.finite(res$weights)))
  expect_equal(sum(res$weights), nrow(df), tolerance = 1e-6)
})

test_that("auto_collapse moves rare mass into __other__ correctly (eb79.17)", {
  df <- data.frame(x = factor(c(rep("a", 100), rep("b", 95), rep("z", 5))))
  res <- harvest(df, list(x = c(a = .45, b = .45, z = .10)),
                 auto_collapse = TRUE, max_weight = 5)

  margins <- tapply(res$weights, df$x, sum) / nrow(df)
  other_share <- margins[["z"]]
  ab_share <- margins[["a"]] + margins[["b"]]

  expect_equal(other_share, .10, tolerance = 1e-2)
  expect_equal(ab_share + other_share, 1, tolerance = 1e-6)
})

test_that("auto_collapse warns when out-of-vocabulary values become NA (eb79.17)", {
  df <- data.frame(x = factor(c(rep("a", 100), rep("b", 95), rep("z", 5), "q")))
  .muffle_sparse_category_diag(expect_warning(
    harvest(df, list(x = c(a = .45, b = .45, z = .10)),
            auto_collapse = TRUE, max_weight = 5),
    "became NA"
  ))
})

test_that("auto_collapse is a no-op and silent when there are no rare/OOV values (eb79.17)", {
  df <- data.frame(x = factor(c(rep("a", 50), rep("b", 50))))
  expect_no_warning(
    res <- harvest(df, list(x = c(a = .5, b = .5)),
                   auto_collapse = TRUE, max_weight = 5)
  )
  expect_equal(attr(res, "result")$status, 0L)   # 0 = RK_OK
})

test_that("collapse_vars alone triggers collapse without auto_collapse (dtkn.2)", {
  # CR-F2: explicit collapse_vars must take effect even though auto_collapse
  # defaults FALSE. The prior gate `!isFALSE(auto_collapse)` made this a silent
  # no-op. "q" is OOV after z collapses into __other__ -> "became NA" warning
  # is emitted ONLY if the collapse path actually ran.
  df <- data.frame(x = factor(c(rep("a", 100), rep("b", 95), rep("z", 5), "q")))
  .muffle_sparse_category_diag(expect_warning(
    harvest(df, list(x = c(a = .45, b = .45, z = .10)),
            collapse_vars = "x", max_weight = 5),
    "became NA"
  ))
})

test_that("collapse_vars naming no target margin warns loudly (dtkn.2)", {
  df <- data.frame(x = factor(c(rep("a", 100), rep("b", 95), rep("z", 5))))
  .muffle_sparse_category_diag(expect_warning(
    harvest(df, list(x = c(a = .45, b = .45, z = .10)),
            collapse_vars = "typo_var", max_weight = 5),
    "no target margin named"
  ))
})
