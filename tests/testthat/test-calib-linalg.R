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
# Gill-Murray-Wright (1981) column-norm bound on LDLT (C1 fix)
# ──────────────────────────────────────────────────────────────────────────────

test_that("calib_linalg: GMW bound prevents L blow-up on nearly-singular margins", {
  # Construct a problem where the unmodified LDLT would form L[i,j] = s/d_j
  # with d_j ~= 0 and s ~= O(1), producing huge L entries (>1e6) that
  # corrupt back-substitution. With the GMW column-norm bound,
  #   d_j >= theta_j^2 / beta^2,
  # so |L[i,j]| <= sqrt(beta^2) <= sqrt(gamma * sqrt(n^2-1)) -- bounded by O(1).
  #
  # We verify the bound indirectly: greg on a nearly-singular system must
  # produce finite weights of moderate magnitude (no >1e6 entries) and
  # bounded margin error. Pre-fix this would either NaN or yield wild weights.
  set.seed(2024)
  n <- 400L
  # Two near-collinear margins: cell counts for category combinations are
  # heavily skewed, making the normal-equation matrix nearly singular.
  data <- data.frame(
    a = factor(sample(c("1", "2"), n, replace = TRUE, prob = c(0.999, 0.001))),
    b = factor(sample(c("1", "2"), n, replace = TRUE, prob = c(0.999, 0.001)))
  )
  target <- list(
    a = c("1" = 0.5, "2" = 0.5),
    b = c("1" = 0.5, "2" = 0.5)
  )

  w <- tryCatch(
    harvest(data, target, method = "greg",
            convergence = list(absolute = 1e-6),
            attach_weights = FALSE),
    error = function(e) e
  )

  # Either greg detects singularity and errors cleanly (C2 guard), or it
  # produces finite, bounded weights (C1 GMW bound). Both are acceptable —
  # the previous behaviour (silent NaN / 1e10-magnitude weights) is not.
  if (inherits(w, "error")) {
    expect_match(conditionMessage(w),
                 "LDLT singular|singular|cell|infeasible|underdetermined",
                 ignore.case = TRUE,
                 info = "near-singular problem must error with a meaningful message, not crash silently")
  } else {
    expect_true(all(is.finite(w)),
                info = "GMW bound: no NaN/Inf weights from near-singular N")
    expect_lt(max(abs(w)), 1e6,
              label = "GMW bound: weights must be bounded (pre-fix would give ~1e10)")
  }
})

# ──────────────────────────────────────────────────────────────────────────────
# CR-A1 (mxcl.1): cholesky_solve must be a TRUE symmetric solve, not a diagonal
# (Jacobi) solve. Regression from commit 61da14b (G2): the normal-equations
# writer fills the row-major LOWER triangle (= column-major UPPER), while dpotrf/
# dpotrs read uplo='L' (= column-major LOWER = zeros off-diagonal) → the solve
# silently degenerates to diag(N)^-1 b.
#
# Discriminator: greg convergence is active-set based (!any_clamped), NOT
# residual based. On an UNBOUNDED, feasible, correlated 2-margin problem one
# full-Newton step (correct X'DX solve) satisfies every margin EXACTLY, so greg
# breaks at iteration 1 with max_error ~ 1e-8. The diagonal-solve bug drops the
# a<->b off-diagonal coupling, so the same iteration-1 break leaves margins
# unmet (max_error >> 1e-3). Both report status=0; max_error separates them.
# ──────────────────────────────────────────────────────────────────────────────

test_that("calib_linalg: cholesky_solve is a full symmetric solve, not diagonal (CR-A1)", {
  set.seed(20260703)
  n <- 600L
  a <- factor(sample(c("1", "2", "3"), n, replace = TRUE, prob = c(0.5, 0.3, 0.2)))
  # b moderately correlated with a: copy a, flip ~20% to a random level so the
  # a<->b off-diagonal blocks of X'DX are clearly non-negligible. (X'DX is exactly
  # rank-deficient by 1 — two complete margins share the sum constraint — but its
  # smallest NONZERO eigenvalue stays large, and the RHS b ⊥ the null vector, so
  # the GMW-regularized correct solve still lands the residual near machine eps.)
  b_chr <- as.character(a)
  flip <- runif(n) < 0.20
  b_chr[flip] <- sample(c("1", "2", "3"), sum(flip), replace = TRUE)
  b <- factor(b_chr, levels = c("1", "2", "3"))
  data <- data.frame(a = a, b = b)

  # Feasible targets shifted a few points off the sample proportions; mild enough
  # that the unbounded Newton weights never touch the loose [0.01, 100] bounds.
  target <- list(
    a = c("1" = 0.45, "2" = 0.33, "3" = 0.22),
    b = c("1" = 0.44, "2" = 0.34, "3" = 0.22)
  )

  w_greg <- harvest(data, target, method = "greg",
                    min_weight = 0.01, max_weight = 100,
                    attach_weights = FALSE)
  r_greg <- attr(w_greg, "result")

  expect_equal(r_greg$status, 0L, info = "greg reports OK (active-set converged)")
  # A TRUE symmetric solve nails every margin in one Newton step; a diagonal
  # (Jacobi) solve drops the a<->b coupling and leaves a large residual here.
  expect_lt(r_greg$max_error, 1e-6)
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
    expect_gt(val[1], 1e-15, label = "schur_nu must be positive (non-degenerate; warm-start from oris may give small but valid schur_nu)")
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
