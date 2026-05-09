test_that("design_effect matches Kish formula", {
  w <- c(1, 2, 3, 4)
  expected <- length(w) * sum(w^2) / sum(w)^2
  expect_equal(design_effect(w), expected, tolerance = 1e-12)
})

test_that("design_effect 4-arg matches Henry-Valliant formula (weighted mean)", {
  w <- c(1, 2, 3, 4)
  y <- c(10, 20, 30, 40)
  # Weighted mean: sum(w*y)/sum(w) = (10+40+90+160)/10 = 30
  y_bar_w <- sum(w * y) / sum(w)
  var_w   <- sum(w * (y - y_bar_w)^2) / sum(w)   # = 100
  var_u   <- var(y)                                 # = 166.667
  expected <- var_w / var_u                         # = 0.6
  expect_equal(design_effect(w, outcome = y), expected, tolerance = 1e-10)
})

test_that("effective_sample_size = n / design_effect", {
  w <- c(1, 2, 3, 4)
  expect_equal(effective_sample_size(w), length(w) / design_effect(w), tolerance = 1e-12)
})
