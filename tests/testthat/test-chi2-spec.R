# test-chi2-spec.R
# Verify spec-correct Pearson chi2: bins with T=0 (or pop<=kMetricEps) must
# contribute ZERO to chi2. Also exercises rare-bin distortion introduced by the
# prior kChi2Floor=1.0 Laplace smoothing (leafblower-p0iy).

test_that("chi2-spec: T=0 bin contributes zero (not delta^2)", {
  set.seed(7)
  n <- 2000L
  # Variable 'a' has 3 categories; category "3" has T=0 (absent from target).
  data <- data.frame(
    a = factor(sample(c("1", "2", "3"), n, replace = TRUE,
                      prob = c(0.4, 0.4, 0.2)))
  )
  # Only categories "1" and "2" are in target; "3" is deliberately absent.
  target <- list(a = c("1" = 0.5, "2" = 0.5))
  w <- leafblower::harvest(data, target, method = "raking",
                            max_iterations = 500,
                            convergence = list(absolute = 1e-8),
                            attach_weights = FALSE)
  result <- attr(w, "result")
  chi2_lb <- result$chi2

  # Spec-correct chi2: computed inline from calibrated weights.
  # Bins with T=0 contribute 0 by Cochran 1954 expected-count rule.
  W   <- sum(w)
  obs <- tapply(w, data$a, sum)  # bucket aggregates
  # Only "1" and "2" are in target; sum-to-1 check
  chi2_spec <- 0
  for (cat in names(target$a)) {
    T_j   <- target$a[[cat]]
    pop_j <- T_j * W
    obs_j <- if (cat %in% names(obs)) obs[[cat]] else 0
    chi2_spec <- chi2_spec + (obs_j - pop_j)^2 / pop_j
  }
  # Category "3" has T=0: its contribution must be 0 in both spec and library.
  # Under the old kChi2Floor=1.0: contribution for cat "3" would be obs["3"]^2 / (0+1).
  obs_3        <- if ("3" %in% names(obs)) obs[["3"]] else 0
  chi2_old_cat3 <- obs_3^2 / (0 + 1.0)  # what kChi2Floor=1.0 injected

  expect_gt(chi2_old_cat3, 0,
    label = "old formula injects positive chi2 for T=0 bin (regression sentinel)")

  expect_equal(chi2_lb, chi2_spec, tolerance = 1e-6,
    label = "library chi2 matches spec-correct Pearson (T=0 bins excluded)")
})

test_that("chi2-spec: rare bin (T~1e-4) not distorted by floor", {
  set.seed(8)
  n <- 5000L
  # Variable 'b' has rare category "2" with T=1e-4 (pop = 1e-4 * n = 0.5).
  # With kChi2Floor=1.0: denom = 0.5+1 = 1.5 → chi2 contribution cut 3x.
  # With spec-correct formula: pop = 0.5, chi2 contribution is correct.
  # prob=5e-3 ensures ~25 obs of "2" in n=5000 — structurally feasible.
  data <- data.frame(
    b = factor(sample(c("1", "2"), n, replace = TRUE,
                      prob = c(1 - 5e-3, 5e-3)))
  )
  T_rare <- 1e-4
  target <- list(b = c("1" = 1 - T_rare, "2" = T_rare))
  w <- leafblower::harvest(data, target, method = "raking",
                            max_iterations = 200,
                            convergence = list(absolute = 1e-8),
                            attach_weights = FALSE)
  result <- attr(w, "result")
  chi2_lb <- result$chi2

  # Spec-correct chi2 inline
  W   <- sum(w)
  obs <- tapply(w, data$b, sum)
  chi2_spec <- 0
  for (cat in names(target$b)) {
    T_j   <- target$b[[cat]]
    pop_j <- T_j * W
    if (pop_j > 1e-10) {
      obs_j <- if (cat %in% names(obs)) obs[[cat]] else 0
      chi2_spec <- chi2_spec + (obs_j - pop_j)^2 / pop_j
    }
  }

  # Regression sentinel: demonstrate the old formula would give a materially
  # different value when obs != pop (e.g., for a hypothetical residual of 1.0
  # unit at the rare bin). Use pop_2 from the converged run.
  # With kChi2Floor=1.0: denom = pop_2 + 1. With spec-correct: denom = pop_2.
  # For pop_2 ~ 0.5: ratio = (0.5+1)/0.5 = 3x difference.
  pop_2      <- T_rare * W  # expected count for rare bin
  obs_2_hyp  <- pop_2 + 1.0    # hypothetical 1-unit residual
  chi2_old_hyp  <- (obs_2_hyp - pop_2)^2 / (pop_2 + 1)
  chi2_spec_hyp <- (obs_2_hyp - pop_2)^2 / pop_2
  expect_false(isTRUE(all.equal(chi2_old_hyp, chi2_spec_hyp, tolerance = 1e-3)),
    label = "old kChi2Floor=1.0 gives different chi2 for rare bin (sentinel)")

  expect_equal(chi2_lb, chi2_spec, tolerance = 1e-6,
    label = "library chi2 matches spec-correct Pearson for rare bin")
})
