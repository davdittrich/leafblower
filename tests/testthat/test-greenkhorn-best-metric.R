library(leafblower)

# CXX.2 (leafblower-5fm8.5): greenkhorn best-iterate must be selected on the
# configured convergence metric (select_metric(cfg.metric, m)), not on the errRp
# fast proxy. For an EXPLICIT metric="max_err" run select_metric == max(errRp), so
# the reported best_error/max_error is byte-identical to the historical code; for
# kl/chi2 (kl is the greenkhorn DEFAULT) the reported value must be on the
# requested metric's scale. Pre-fix the slot held the errRp magnitude mislabeled
# as the configured metric — visible on a NON-converging kl run (see below).

make_feasible_df <- function() {
  set.seed(7)
  n <- 200L
  data.frame(
    a = factor(sample(c("1", "2", "3"), n, TRUE)),
    b = factor(sample(c("x", "y"), n, TRUE))
  )
}
feasible_targets <- function() {
  list(a = c("1" = 1 / 3, "2" = 1 / 3, "3" = 1 / 3),
       b = c(x = 0.5, y = 0.5))
}

test_that("CXX.2: greenkhorn default-metric feasible problem converges (no regression)", {
  # Default convergence metric for greenkhorn is "kl"; on a fully feasible problem
  # KL -> 0 as it converges, so best_error/max_error collapse to ~0 regardless of
  # scale. This is a no-regression smoke test (status 0, tiny residual).
  w <- harvest(make_feasible_df(), feasible_targets(),
    method = "greenkhorn", max_iterations = 500L, attach_weights = FALSE)
  r <- attr(w, "result")
  expect_equal(r$status, 0L)
  expect_lt(r$max_error, 1e-3)
  expect_equal(r$max_error, r$best_error, tolerance = 1e-12)
})

test_that("CXX.2: greenkhorn metric=kl reports best_error on the KL scale", {
  w <- harvest(make_feasible_df(), feasible_targets(),
    method = "greenkhorn",
    convergence = list(metric = "kl", rule = "threshold", tol = 1e-8),
    max_iterations = 500L, attach_weights = FALSE)
  r <- attr(w, "result")
  # best-iterate selected via select_metric(KL, m): max_error == best_error and
  # both are the KL value of the best iterate (finite, non-negative). Pre-fix
  # this slot held the errRp max instead of the KL value.
  expect_true(is.finite(r$best_error))
  expect_gte(r$best_error, 0)
  expect_equal(r$max_error, r$best_error, tolerance = 1e-12)
  # The convergence metric recorded must be KL, confirming the configured
  # metric drove the run that the best-iterate is now selected against.
  expect_equal(r$convergence_used$metric, "kl")
})

# Recompute the internal KL metric (compute_cell_metrics in calib_dispatch.hpp)
# EXACTLY: m.kl = max_k Σ_j T_kj · log((T_kj + eps)/(S_p_kj + eps)), eps=1e-10,
# S_p_kj = (Σ weights in bucket k,j) / Σ weights. This is the scale select_metric
# returns for metric="kl"; the returned greenkhorn weights are the best-iterate
# snapshot, so recomputing on them must reproduce r$best_error.
kKlMetricEps <- 1e-10
internal_kl_max <- function(weights, df, tgt) {
  W <- sum(weights)
  kl_max <- 0
  for (nm in names(tgt)) {
    lv <- names(tgt[[nm]])
    S <- tapply(weights, df[[nm]], sum)[lv]
    S[is.na(S)] <- 0
    tk <- unname(tgt[[nm]])
    Sr <- as.numeric(S / W)
    kl_k <- sum(ifelse(tk > 0, tk * log((tk + kKlMetricEps) / (Sr + kKlMetricEps)), 0))
    kl_max <- max(kl_max, kl_k)
  }
  kl_max
}
errRp_max <- function(weights, df, tgt) {
  W <- sum(weights)
  e <- 0
  for (nm in names(tgt)) {
    lv <- names(tgt[[nm]])
    S <- tapply(weights, df[[nm]], sum)[lv]
    S[is.na(S)] <- 0
    tk <- unname(tgt[[nm]])
    Sr <- as.numeric(S / W)
    e <- max(e, max(abs(Sr - tk)))
  }
  e
}

test_that("CXX.2: NON-converging metric=kl reports best_error on KL scale, not errRp proxy", {
  # Skewed targets + tiny iteration budget => greenkhorn cannot converge (status
  # != 0). This is where the scale-mix bug manifested: pre-fix, best.best_metric
  # was populated per-iteration on the errRp (max marginal-residual) scale, so the
  # reported best_error held an errRp value mislabeled as KL. Post-fix the reported
  # best-iterate lives strictly on select_metric(cfg.metric=KL, m).
  set.seed(7)
  n <- 200L
  df <- data.frame(
    a = factor(sample(c("1", "2", "3"), n, TRUE)),
    b = factor(sample(c("x", "y"), n, TRUE))
  )
  tgt <- list(a = c("1" = 0.6, "2" = 0.3, "3" = 0.1), b = c(x = 0.7, y = 0.3))

  w <- suppressWarnings(harvest(df, tgt,
    method = "greenkhorn",
    convergence = list(metric = "kl", rule = "threshold", tol = 1e-12),
    max_iterations = 4L, attach_weights = FALSE))
  r <- attr(w, "result")
  ww <- as.numeric(w)

  # Did NOT converge — the regime where the bug is observable.
  expect_false(identical(r$status, 0L))

  kl_ref    <- internal_kl_max(ww, df, tgt)
  errRp_ref <- errRp_max(ww, df, tgt)

  # The two scales must be genuinely DIFFERENT here, otherwise the assertion
  # below could not distinguish the fix from the bug (degenerate guard).
  expect_gt(abs(kl_ref - errRp_ref), 1e-3)

  # best_error / max_error live on the KL scale: they equal the internal KL of
  # the returned (best-iterate) weights to machine precision, and do NOT equal
  # the errRp marginal-residual magnitude (which is what the pre-fix code stored).
  expect_equal(r$best_error, kl_ref, tolerance = 1e-12)
  expect_equal(r$max_error, r$best_error, tolerance = 1e-12)
  expect_gt(abs(r$best_error - errRp_ref), 1e-3)
})
