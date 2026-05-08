test_that("design_effect matches Kish formula", {
  w <- c(1, 2, 3, 4)
  expected <- length(w) * sum(w^2) / sum(w)^2
  expect_equal(design_effect(w), expected, tolerance = 1e-12)
})

test_that("design_effect 4-arg matches Henry-Valliant calibration design effect (Cochran 1977 weighted var)", {
  w <- c(1, 2, 3, 4)
  y <- c(10, 20, 30, 40)
  y_bar_w <- sum(w * y) / sum(w)
  # Cochran (1977) §4.5 unbiased weighted variance — reduces to Bessel for uniform w.
  var_w_denom <- (sum(w)^2 - sum(w^2)) / sum(w)
  var_w   <- sum(w * (y - y_bar_w)^2) / var_w_denom
  var_u   <- var(y)
  expected <- var_w / var_u
  expect_equal(design_effect(w, outcome = y), expected, tolerance = 1e-10)
})

test_that("design_effect 4-arg returns 1.0 for uniform weights (Cochran reduction)", {
  w <- c(1, 1, 1, 1, 1)
  y <- c(10, 20, 30, 40, 50)
  expect_equal(design_effect(w, outcome = y), 1.0, tolerance = 1e-12)
})

test_that("effective_sample_size = n / design_effect", {
  w <- c(1, 2, 3, 4)
  expect_equal(effective_sample_size(w), length(w) / design_effect(w), tolerance = 1e-12)
})
