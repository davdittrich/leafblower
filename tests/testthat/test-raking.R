test_that("iEPPA converges: 1 margin, 2 cats, no bounds", {
  set.seed(42)
  n   <- 100L
  df  <- data.frame(x = factor(sample(c("a","b"), n, replace=TRUE, prob=c(0.7,0.3))))
  tgt <- list(x = c(a=0.5, b=0.5))
  result <- leafblower::harvest(df, tgt, method = "raking", convergence = list(absolute = 1e-6))
  expect_true(attr(result, "algorithm") == "raking")
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("iEPPA respects max_weight=2 on tight bounds", {
  set.seed(7)
  n   <- 10000L
  df  <- data.frame(
    age = factor(sample(c("Y","M","O"), n, replace=TRUE, prob=c(0.40,0.33,0.27))),
    sex = factor(sample(c("M","F"), n, replace=TRUE, prob=c(0.55,0.45))),
    edu = factor(sample(c("HS","Col","Grad"), n, replace=TRUE, prob=c(0.4,0.4,0.2)))
  )
  tgt <- list(
    age = c(Y=0.33, M=0.34, O=0.33),
    sex = c(M=0.50, F=0.50),
    edu = c(HS=0.35, Col=0.45, Grad=0.20)
  )
  result <- leafblower::harvest(df, tgt, method = "raking", max_weight=2, convergence = list(absolute = 1e-6))
  expect_true(max(result$weights) <= 2.0 + 1e-10)
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("iEPPA respects min_weight=0.5", {
  set.seed(3)
  n   <- 10000L
  df  <- data.frame(
    x = factor(sample(c("a","b","c","d","e"), n, replace=TRUE))
  )
  tgt <- list(x = c(a=0.2, b=0.2, c=0.2, d=0.2, e=0.2))
  result <- leafblower::harvest(df, tgt, method = "raking", min_weight=0.5, max_weight=5, convergence = list(absolute = 1e-6))
  expect_true(min(result$weights) >= 0.5 - 1e-10)
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})

test_that("iEPPA output weights have mean=1 and respect bounds", {
  set.seed(5L)
  n   <- 1000L
  df  <- data.frame(x = factor(sample(c("a", "b", "c"), n, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.3, c = 0.2))
  res <- leafblower::harvest(df, tgt, method = "raking",
                              max_weight = 2.0, min_weight = 0.2,
                              convergence = list(absolute = 1e-6),
                              attach_weights = FALSE)
  # mean=1 is guaranteed by both the old fixup loop and the new Dykstra projection;
  # this test guards against regressions in the P2 refactor.
  expect_equal(mean(res), 1.0, tolerance = 1e-10)
  expect_true(max(res) <= 2.0 + 1e-10)
  expect_true(min(res) >= 0.2 - 1e-10)
})

test_that("descent monitor aborts early on stalled errRp trajectory", {
  # Input empirically validated to stall: n=1000, 95/5 class split, max_weight=1.2.
  # Prior measurement (leafblower-370 probe session 2026-04-23): raking hit
  # max_iter=500 with no convergence and 500 consecutive curvature rejections.
  # Monitor fires at 5 consecutive stalled error-checks (50 iters) — well
  # before iter 500.
  set.seed(91)
  n <- 1000
  df <- data.frame(cat = sample(c("A", "B"), n, replace = TRUE, prob = c(0.05, 0.95)))
  tgt <- list(cat = c(A = 0.95, B = 0.05))
  # Emit verbose=1 so we can grep for the monitor message.
  # CalibState.log() routes through Rprintf (src/r_bridge.cpp:21) which writes
  # to stdout, not R's message sink. Use capture.output(type = "output").
  t0 <- Sys.time()
  msgs <- capture.output(
    tryCatch(
      suppressWarnings(leafblower::harvest(df, tgt, method = "raking",
                                max_weight = 1.2,
                                max_iterations = 500,
                                verbose = 1L,
                                convergence = list(absolute = 1e-6))),
      error = function(e) NULL  # B9: infeasible error expected; msgs still captured
    ),
    type = "output"
  )
  elapsed <- as.numeric(Sys.time() - t0, units = "secs")
  # Wall-clock guard: monitor should trigger early; whole call must be under 5s
  # even if n_max_iter=500 is reached without the monitor.
  expect_lt(elapsed, 5)
  # Monitor message matches the exact phrase emitted by src/raking.cpp.
  probe <- paste(msgs, collapse = "\n")
  expect_match(probe, "stalled for [0-9]+ consecutive checks",
               info = paste("expected descent-monitor message; got:", probe))
})

test_that("raking with binding capacity constraints returns best achievable result (not INFEAS)", {
  # n=20, 2 obs in cat "a", 18 in cat "b".
  # Target: 90% in "a" → requires average weight 9 per obs in "a".
  # max_weight=1.5 caps this: cell mass for "a" is at most 1.5*2=3 << 18.
  # This is a binding capacity constraint, NOT structural infeasibility.
  # Correct behavior: return the best achievable constrained optimum (STALL or BUDGET),
  # not INFEAS. INFEAS is reserved for lower-bound violations (min_weight sum > target).
  result <- suppressWarnings(
    leafblower::harvest(data.frame(x = factor(c(rep("a", 2L), rep("b", 18L)))),
            target = list(x = c(a=0.9, b=0.1)),
            method = "raking",
            min_weight = 0.5, max_weight = 1.5,
            max_iterations = 200L,
            convergence = list(absolute = 1e-4),
            attach_weights = FALSE)
  )
  s <- attr(result, "result")$status
  # Must NOT be INFEAS (2) — capacity capping is normal bounded calibration.
  expect_false(s == 2L, info = paste("raking returned INFEAS on a capacity-constrained (not infeasible) problem; status:", s))
  # Must be STALL (5) or BUDGET (4) — best achievable constrained result.
  expect_true(s %in% c(4L, 5L), info = paste("expected STALL(5) or BUDGET(4), got:", s))
})

test_that("B16: SQUAREM stall detection does not fire spuriously after fallback", {
  result <- leafblower::harvest(
    data.frame(x = rep(c("A","B"), 50), w=1),
    target = list(x = c(A=0.5, B=0.5)),
    method = "raking",
    convergence = list(rule="improvement", pct=1e-4),
    accelerate = TRUE,
    max_iterations = 200L
  )
  # A uniform problem with target=empirical proportion converges quickly.
  expect_lt(attr(result,"result")$iterations, 20L)
})

test_that("R8: accelerate=TRUE with greedy scheduler runs without error", {
  expect_no_error({
    result <- leafblower::harvest(
      data.frame(x = rep(c("A","B"), 25)),
      target = list(x = c(A=0.5, B=0.5)),
      method = "raking",
      convergence = list(rule="improvement", pct=1e-4),
      accelerate = TRUE,
      scheduler  = "greedy",
      max_iterations = 50L
    )
    expect_gt(attr(result,"result")$iterations, 0L)
  })
})
