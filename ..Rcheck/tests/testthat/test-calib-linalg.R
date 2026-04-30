# Direct tests for LDLT factorization and compute_normal_equations
# Tests the implementation indirectly via greg solver (which uses LDLT internally).
# No direct C-level test bridge; these verify LDLT correctness by checking
# that greg achieves low margin errors on synthetic balanced problems.

test_that("calib_linalg: LDLT correctness via greg on balanced 1-margin problem", {
  # Simple balanced problem: 3-level factor, uniform distribution, balanced target
  # greg must solve the normal equations (computed via compute_normal_equations)
  # using LDLT factorization. If LDLT is wrong, convergence fails or error is large.
  set.seed(42)
  n <- 300L
  data <- data.frame(a = factor(sample(c("1", "2", "3"), n, TRUE)))
  target <- list(a = c("1" = 0.33, "2" = 0.34, "3" = 0.33))

  w_greg <- harvest(data, target, method = "greg",
                    convergence = list(absolute = 1e-6),
                    attach_weights = FALSE)
  r_greg <- attr(w_greg, "result")

  # greg should converge (status 0) and achieve chi2 < 0.01
  expect_equal(r_greg$status, 0L, info = "greg must converge")
  expect_true(r_greg$chi2 < 0.01, info = "greg chi2 must be small")
  # Margin errors should be near zero (LDLT solved correctly)
  expect_true(r_greg$max_error < 0.01, info = "greg max margin error < 1%")
})

test_that("calib_linalg: LDLT handles near-singular (zero category) via regularization", {
  # Skewed distribution: category "3" appears only ~2% of the time
  # The normal matrix N may have a near-zero diagonal for "3"
  # Gill-Murray perturbation (if implemented) or LDLT robustness must prevent division by zero
  set.seed(99)
  n <- 500L
  data <- data.frame(
    a = factor(sample(c("1", "2", "3"), n, replace = TRUE,
                      prob = c(0.49, 0.49, 0.02)))
  )
  target <- list(a = c("1" = 0.45, "2" = 0.45, "3" = 0.1))

  w_greg <- harvest(data, target, method = "greg",
                    convergence = list(absolute = 1e-6),
                    attach_weights = FALSE)
  r_greg <- attr(w_greg, "result")

  # Should not crash or return Inf/NaN
  expect_true(is.finite(r_greg$max_error),
              info = "no Inf/NaN from near-singular N")
  expect_true(r_greg$status %in% c(0L, 4L, 5L),
              info = "status 0 (converged) or 4 (BUDGET) or 5 (STALL), not crash")
})

test_that("calib_linalg: LDLT on 2-margin balanced problem", {
  # Two-way problem: age × sex
  # Normal equations are larger (3×3 + 2×2 - 1 variables), tests LDLT at scale
  set.seed(123)
  n <- 1000L
  data <- data.frame(
    age = factor(sample(c("Y", "M", "O"), n, replace = TRUE,
                        prob = c(0.60, 0.30, 0.10))),
    sex = factor(sample(c("M", "F"), n, replace = TRUE,
                        prob = c(0.70, 0.30)))
  )
  target <- list(
    age = c(Y = 0.33, M = 0.40, O = 0.27),
    sex = c(M = 0.49, F = 0.51)
  )

  w_greg <- harvest(data, target, method = "greg",
                    convergence = list(absolute = 1e-6),
                    attach_weights = FALSE)
  r_greg <- attr(w_greg, "result")

  # Normal equations are 4×4 (5 parameters, 1 sum constraint)
  # LDLT factorization must work correctly
  expect_equal(r_greg$status, 0L, info = "greg converges on 2-margin problem")
  expect_true(r_greg$max_error < 0.01, info = "max margin error < 1%")
})

# ──────────────────────────────────────────────────────────────────────────────
# chebyshev nu-fix: reference elimination makes schur_nu > 0
# ──────────────────────────────────────────────────────────────────────────────

test_that("chebyshev: schur_nu diagnostic logged at verbose=2 (non-degeneracy)", {
  set.seed(42); n <- 200L
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, TRUE)),
    b = factor(sample(c("1","2"),     n, TRUE))
  )
  target <- list(a = c("1"=0.4,"2"=0.4,"3"=0.2), b = c("1"=0.6,"2"=0.4))

  # Capture Rprintf output (stdout) via capture.output(type="output")
  log_lines <- capture.output(
    leafblower::harvest(data, target, method = "chebyshev",
                        min_weight = 0.2, max_weight = 5,
                        max_iterations = 5, attach_weights = FALSE,
                        verbose = 2),
    type = "output"
  )
  schur_lines <- grep("schur_nu=", log_lines, value = TRUE)

  expect_gt(length(schur_lines), 0L,
            label = "schur_nu should be logged at verbose=2")
  if (length(schur_lines) > 0) {
    val <- as.numeric(regmatches(schur_lines[1],
                                 regexpr("[0-9]+\\.?[0-9]*[eE][+-]?[0-9]+|[0-9]+\\.?[0-9]*",
                                         schur_lines[1])))
    expect_gt(val[1], 1e-6, label = "schur_nu must be positive (non-degenerate)")
  }
})

test_that("chebyshev: single-category margin does not crash", {
  set.seed(99); n <- 100L
  data <- data.frame(
    a = factor(sample(c("1","2"), n, TRUE)),
    b = factor(rep("x", n))   # single category — trivial constraint
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("x"=1.0))
  r <- leafblower::harvest(data, target, method = "chebyshev",
                           min_weight = 0.2, max_weight = 5,
                           max_iterations = 20, attach_weights = FALSE,
                           verbose = 0)
  expect_true(all(is.finite(r)),
              label = "chebyshev with 1-cat margin: weights must be finite")
  expect_true(all(r >= 0),
              label = "chebyshev with 1-cat margin: weights must be non-negative")
})
