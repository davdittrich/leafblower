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
  expect_warning(
    harvest(df, list(x = c(a = .45, b = .45, z = .10)),
            auto_collapse = TRUE, max_weight = 5),
    "became NA"
  )
})

test_that("auto_collapse is a no-op and silent when there are no rare/OOV values (eb79.17)", {
  df <- data.frame(x = factor(c(rep("a", 50), rep("b", 50))))
  expect_no_warning(
    res <- harvest(df, list(x = c(a = .5, b = .5)),
                   auto_collapse = TRUE, max_weight = 5)
  )
  expect_equal(attr(res, "result")$status, 0L)   # 0 = RK_OK
})
