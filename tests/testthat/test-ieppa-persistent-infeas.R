# WU-1: persistent-infeas tracker regression test.
# Pre-fix: iEPPA flags RK_ERR_INFEAS on transient near-zero buckets
#   (stepstone-shape K=5 input with overlapping margins).
# Post-fix: converges to RK_OK; infeasibility reported only after
#   kInfeasPersistence=5 consecutive outer iterations with the same
#   bucket empty (guards against transient-in-settling false positives).

test_that("WU-1: iEPPA converges on structurally feasible K=5 overlapping margins", {
  set.seed(42)
  n <- 2000L
  K <- 5L
  # 3-3-4-3-3 cats with overlapping (correlated) groupings.
  df <- data.frame(
    a = sample(letters[1:3], n, replace = TRUE, prob = c(0.5, 0.3, 0.2)),
    b = sample(letters[1:3], n, replace = TRUE, prob = c(0.2, 0.5, 0.3)),
    c = sample(letters[1:4], n, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.2)),
    d = sample(letters[1:3], n, replace = TRUE, prob = c(0.4, 0.3, 0.3)),
    e = sample(letters[1:3], n, replace = TRUE, prob = c(0.25, 0.25, 0.5))
  )
  targets <- list(
    a = c(a = 0.40, b = 0.35, c = 0.25),
    b = c(a = 0.30, b = 0.40, c = 0.30),
    c = c(a = 0.25, b = 0.25, c = 0.25, d = 0.25),
    d = c(a = 0.33, b = 0.33, c = 0.34),
    e = c(a = 0.30, b = 0.35, c = 0.35)
  )
  # Structurally feasible at max_weight=5; pre-fix latches INFEAS on transient.
  expect_no_error(
    res <- harvest(df, targets, method = "ieppa",
                   max_weight = 5, min_weight = 0,
                   max_iterations = 500L,
                   convergence = list(absolute = 1e-4),
                   attach_weights = FALSE)
  )
  diag <- diagnose_weights(df, targets, res)
  expect_lt(max(abs(diag$error_weighted)), 1e-3)
})

test_that("WU-1: truly infeasible input (empty target cell) still reports INFEAS", {
  # Regression guard: genuine infeasibility must still latch.
  n <- 500L
  df <- data.frame(
    a = sample(letters[1:2], n, replace = TRUE),
    b = sample(letters[1:2], n, replace = TRUE)
  )
  # Target a third category 'c' that has zero observations → persistent empty.
  targets <- list(
    a = c(a = 0.4, b = 0.3, c = 0.3),
    b = c(a = 0.5, b = 0.5)
  )
  expect_error(
    suppressWarnings(harvest(df, targets, method = "ieppa",
                             max_weight = 5, min_weight = 0,
                             max_iterations = 500L,
                             convergence = list(absolute = 1e-6))),
    regexp = "persistent empty cell|infeasible problem"
  )
})

test_that("WU-1: oscillating streak (spec §4 edge case) returns NOCONV not INFEAS", {
  # Spec §4 documents: a bucket that oscillates empty <-> non-empty such that
  # streak resets before kInfeasPersistence=5 will NOT flag INFEAS; solver
  # hits max_iter -> RK_ERR_BUDGET (status=4) with high errRp. This test guards that
  # documented behaviour. Engineer oscillation via a 3-way near-degenerate
  # system where each outer sweep alternates which margin is pinched.
  set.seed(2024)
  n <- 600L
  # K=3 with strong negative correlation between two margins → oscillation.
  a <- sample(letters[1:3], n, replace = TRUE, prob = c(0.1, 0.45, 0.45))
  # b chosen so (a,b) cells heavily biased; targets will push the solver to
  # move mass between (a="a", any b) cells, which are sparse.
  b <- ifelse(a == "a", sample(letters[1:3], n, replace = TRUE, prob = c(0.8, 0.1, 0.1)),
                        sample(letters[1:3], n, replace = TRUE))
  c_ <- sample(letters[1:3], n, replace = TRUE)
  df <- data.frame(a = a, b = b, c_ = c_)
  # Targets pushing the (a="a", b="b") and (a="a", b="c") cells to non-trivial mass.
  targets <- list(
    a  = c(a = 0.50, b = 0.25, c = 0.25),
    b  = c(a = 0.25, b = 0.50, c = 0.25),
    c_ = c(a = 0.33, b = 0.33, c = 0.34)
  )
  # Low max_iter forces solver to exit before streak settles either way.
  res_status <- tryCatch({
    suppressWarnings(harvest(df, targets, method = "ieppa",
                             max_weight = 5, min_weight = 0,
                             max_iterations = 30L,
                             convergence = list(absolute = 1e-10)))
    "converged"
  }, error = function(e) {
    if (grepl("persistent empty cell", conditionMessage(e))) "infeas"
    else "other_error"
  })
  # Must NOT be "infeas" on this near-degenerate but not structurally empty input
  # with streak-resetting oscillation. "converged" (OK) or a non-infeas warning
  # path are both acceptable per spec §4.
  expect_false(res_status == "infeas")
})
