test_that("F(0) = 1 for logit link", {
  # Calls C function that doesn't exist yet — must FAIL
  expect_equal(.Call("C_logit_F_at_zero", 0.5, 5.0), 1.0, tolerance = 1e-12)
})

test_that("F(u) stays in [L,U] for logit link", {
  # Also RED until implemented
  vals <- .Call("C_logit_range_check", 0.5, 5.0, as.double(seq(-10, 10, by=0.5)))
  expect_true(all(vals >= 0.5 - 1e-12))
  expect_true(all(vals <= 5.0 + 1e-12))
})

test_that("H prime equals F for logit link (numerical diff)", {
  # RED
  result <- .Call("C_logit_Hprime_check", 0.5, 5.0, 1.0)
  expect_equal(result, 0.0, tolerance = 1e-8)
})

test_that("exp link: F(u) = exp(u)", {
  # RED
  expect_equal(.Call("C_logit_F_at_zero", 0.0, Inf), 1.0, tolerance = 1e-12)
})
