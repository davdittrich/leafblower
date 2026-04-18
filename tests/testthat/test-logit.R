test_that("F(0) = 1 for logit link", {
  expect_equal(.Call("C_logit_F_at_zero", 0.5, 5.0), 1.0, tolerance = 1e-12)
})

test_that("F(u) stays in [L,U] for logit link", {
  vals <- .Call("C_logit_range_check", 0.5, 5.0, as.double(seq(-10, 10, by=0.5)))
  expect_true(all(vals >= 0.5 - 1e-12))
  expect_true(all(vals <= 5.0 + 1e-12))
})

test_that("H prime equals F for logit link (numerical diff)", {
  result <- .Call("C_logit_Hprime_check", 0.5, 5.0, 1.0)
  expect_equal(result, 0.0, tolerance = 1e-8)
})

test_that("exp link: F(u) = exp(u)", {
  expect_equal(.Call("C_logit_F_at_zero", 0.0, Inf), 1.0, tolerance = 1e-12)
})

test_that("H'(u) = F(u) holds near safe_exp clamp boundary (u=559 for L=0,U=5)", {
  # logit_scale = (5-0)/((5-1)*(1-0)) = 1.25
  # logit_scale * 559 = 698.75 (just below the safe_exp clamp at 700)
  # safe_exp clamp preserves H'(u) = F(u) because both F and H use the same
  # safe_exp(logit_scale*u) value; the algebraic identity is maintained.
  result <- .Call("C_logit_Hprime_check", 0.0, 5.0, 559.0)
  expect_equal(result, 0.0, tolerance = 1e-4)  # wider tol: finite-diff truncation error near saturation (F(u)≈U=5)
})
