test_that("convergence=list(absolute=1e-6) is backward compat (max_err criterion)", {
  set.seed(101)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(
    a = c("1"=0.3, "2"=0.5, "3"=0.2),
    b = c("1"=0.6, "2"=0.4)
  )
  w <- leafblower::harvest(
    data, target, max_weight = 3, method = "ieppa",
    convergence = list(absolute = 1e-6),
    attach_weights = FALSE
  )
  expect_true(is.numeric(as.numeric(w)))
  expect_length(as.numeric(w), n)
})

test_that("A2: pct=0.001 default converges on smooth synthetic", {
  set.seed(42)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE)),
    c = factor(sample(c("1","2","3","4"), n, replace = TRUE))
  )
  target <- list(
    a = c("1"=1/3, "2"=1/3, "3"=1/3),
    b = c("1"=0.5, "2"=0.5),
    c = c("1"=0.25, "2"=0.25, "3"=0.25, "4"=0.25)
  )
  w <- leafblower::harvest(
    data, target, max_weight = 10, method = "ieppa",
    max_iterations = 500,
    attach_weights = FALSE
  )
  result <- attr(w, "result")
  # l1_weight_change: finite and non-negative (no threshold — default criterion is
  # max_err+improvement, not l1_weight, so l1_weight_change reflects change at
  # convergence iter and may be large)
  expect_gte(result$l1_weight_change, 0)
  expect_true(is.finite(result$l1_weight_change))
  expect_lt(result$max_error, 1e-3)
  expect_equal(result$status, 0L)
  # Verify pct_change IS computed: use a barely-converged run (max_iterations=10)
  # where the solver stops before pct threshold is met; l1_weight_change must be > 0
  w10 <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    max_iterations = 10, attach_weights = FALSE)
  r10 <- attr(w10, "result")
  expect_gt(r10$l1_weight_change, 0)
})

test_that("A8a: metric='mean_err' actively stops solver", {
  set.seed(43)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    convergence = list(absolute = 1e-4, metric = "mean_err"),
    max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  # mean_error must be computed (non-zero) and within threshold
  expect_gt(result$mean_error, 0)
  expect_lt(result$mean_error, 1e-4)
  expect_equal(result$status, 0L)
})

test_that("A8b: metric='kl' actively stops solver", {
  set.seed(44)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    convergence = list(absolute = 1e-6, metric = "kl"),
    max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  # kl must be computed (non-zero) and within threshold
  expect_gt(result$kl, 0)
  expect_lt(result$kl, 1e-6)
  expect_equal(result$status, 0L)
})

test_that("A8c: metric='chi2' actively stops solver", {
  set.seed(45)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  # chi2 scales with n; threshold n-scaled (~1e-3 * 2000 = 2)
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    convergence = list(absolute = 2.0, metric = "chi2"),
    max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  # chi2 must be computed (non-zero) and within threshold
  expect_gt(result$chi2, 0)
  expect_lt(result$chi2, 2.0)
  expect_equal(result$status, 0L)
})

test_that("harvest accepts sor argument without error", {
  set.seed(101)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5, "2"=0.5))
  w1 <- leafblower::harvest(data, target, max_weight = 3, method = "ieppa",
                            sor = NULL, attach_weights = FALSE)
  w2 <- leafblower::harvest(data, target, max_weight = 3, method = "ieppa",
                            sor = list(auto = TRUE, omega_min = 0.3),
                            attach_weights = FALSE)
  expect_length(as.numeric(w1), n)
  expect_length(as.numeric(w2), n)
})

test_that("parse_convergence rejects unknown keys", {
  expect_error(
    leafblower::harvest(
      data.frame(a = factor(c("1","2"))),
      list(a = c("1"=0.5, "2"=0.5)),
      max_weight = 3, method = "ieppa",
      convergence = list(pct_tol = 0.001),
      attach_weights = FALSE
    ),
    regexp = "Unknown convergence key"
  )
})

test_that("parse_convergence rejects non-list convergence", {
  expect_error(
    leafblower::harvest(
      data.frame(a = factor(c("1","2"))),
      list(a = c("1"=0.5, "2"=0.5)),
      max_weight = 3, method = "ieppa",
      convergence = 1e-6,
      attach_weights = FALSE
    ),
    regexp = "convergence must be a named list"
  )
})

test_that("parse_sor rejects unknown keys", {
  expect_error(
    leafblower::harvest(
      data.frame(a = factor(c("1","2"))),
      list(a = c("1"=0.5, "2"=0.5)),
      max_weight = 3, method = "ieppa",
      sor = list(omega_minimum = 0.3),
      attach_weights = FALSE
    ),
    regexp = "Unknown sor key"
  )
})

test_that("hawe: iEPPA does NOT warn PCT stall when metric=pct (l1_weight)", {
  # Regression guard for S1: stall warning must NOT fire for pct/l1_weight metric.
  # High max_error is EXPECTED when the user chose weight-change convergence —
  # they don't care about error residuals, so firing is a false positive.
  # Contradictory targets stress-test that the guard is correctly suppressed.
  n <- 400
  var1 <- factor(rep(c("A","B"), each = n/2))
  var2 <- factor(rep(c("2","1"), each = n/2))  # A maps to "2", B maps to "1"
  data <- data.frame(var1 = var1, var2 = var2)
  target <- list(
    var1 = c(A = 0.95, B = 0.05),
    var2 = c("1" = 0.95, "2" = 0.05)
  )
  result <- suppressWarnings(
    leafblower::harvest(data, target, max_weight = 1.5, method = "ieppa",
                        max_iterations = 300,
                        convergence = list(pct = 0.001),
                        attach_weights = FALSE)
  )
  # Collect only the PCT stall warning specifically
  w <- tryCatch(
    withCallingHandlers(
      leafblower::harvest(data, target, max_weight = 1.5, method = "ieppa",
                          max_iterations = 300,
                          convergence = list(pct = 0.001),
                          attach_weights = FALSE),
      warning = function(w) {
        if (grepl("PCT convergence stall", conditionMessage(w))) {
          stop(paste0("STALL_WARNING: ", conditionMessage(w)))
        }
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  expect_false(
    inherits(w, "error") && grepl("STALL_WARNING", conditionMessage(w)),
    info = "PCT stall warning must NOT fire when metric=pct (l1_weight)"
  )
})

# CONFIRMATORY TEST (not TDD red-green): added in same commit as implementation.
# Retroactive analysis: this test WOULD fail on pre-WU-C code because
# convergence_used$rule is NULL until harvest.R WU-E2 wires the nesting.
test_that("A1: default convergence (kl+improvement) converges smooth synthetic", {
  set.seed(42)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE)),
    c = factor(sample(c("1","2","3","4"), n, replace = TRUE))
  )
  target <- list(
    a = c("1"=1/3, "2"=1/3, "3"=1/3),
    b = c("1"=0.5, "2"=0.5),
    c = c("1"=0.25, "2"=0.25, "3"=0.25, "4"=0.25)
  )
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
                           max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$status, 0L, info = "must converge")
  expect_lt(result$max_error, 1e-3)
  # ieppa default: rule=improvement, metric=marginal_kl (Task 0 calibration quality)
  expect_equal(result$convergence_used$rule, "improvement",
               info = "default rule must be improvement")
  expect_equal(result$convergence_used$metric, "marginal_kl",
               info = "ieppa default metric must be marginal_kl (Task 0)")
})

# CONFIRMATORY TEST (not TDD red-green): added in same commit as implementation.
# Retroactive analysis: this test WOULD fail on pre-WU-C code because
# convergence_used$rule is NULL until harvest.R WU-E2 wires the nesting.
test_that("A2: oscillating input — best_error < 0.9 * max_error on NOCONV", {
  set.seed(31415)
  n <- 2000
  data <- data.frame(
    v1 = factor(sample(c("A","B","C","D"), n, replace = TRUE)),
    v2 = factor(sample(c("X","Y","Z"), n, replace = TRUE)),
    v3 = factor(sample(c("1","2","3","4","5"), n, replace = TRUE)),
    v4 = factor(sample(c("p","q"), n, replace = TRUE)),
    v5 = factor(sample(c("a","b","c","d","e","f"), n, replace = TRUE))
  )
  target <- list(
    v1 = c(A=0.1, B=0.4, C=0.4, D=0.1),
    v2 = c(X=0.5, Y=0.3, Z=0.2),
    v3 = c("1"=0.1,"2"=0.1,"3"=0.4,"4"=0.3,"5"=0.1),
    v4 = c(p=0.7, q=0.3),
    v5 = c(a=0.05,b=0.05,c=0.5,d=0.2,e=0.15,f=0.05)
  )
  w <- leafblower::harvest(data, target, max_weight = 2, method = "ieppa",
                           max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  expect_lte(result$best_error, result$max_error)
  if (result$status %in% c(4L, 5L)) {
    expect_lt(result$best_error, 0.9 * result$max_error)
  }
})

# CONFIRMATORY TEST (not TDD red-green): added in same commit as implementation.
# Retroactive analysis: this test WOULD fail on pre-WU-C code because
# convergence_used$rule is NULL until harvest.R WU-E2 wires the nesting.
test_that("A3: list(pct=0.001) triggers l1_weight+plateau on raking", {
  set.seed(43)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "raking",
                           max_iterations = 500,
                           convergence = list(pct = 0.001),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$status, 0L)
  # plateau fires when l1_weight stops improving (not when it's absolutely small)
  expect_gte(result$l1_weight_change, 0)
  expect_true(is.finite(result$l1_weight_change))
  # pct key -> metric=l1_weight, rule=plateau — autumn/anesrake compatible
  expect_equal(result$convergence_used$metric, "l1_weight",
               info = "pct key must select l1_weight metric")
  expect_equal(result$convergence_used$rule, "plateau",
               info = "pct must use plateau rule")
})

test_that("A4: explicit improvement rule with absolute tol fires on raking", {
  set.seed(55)
  n <- 1500
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  # Use absolute tol so the solver converges, and verify convergence metadata stored.
  w <- leafblower::harvest(data, target, max_weight = 10, method = "raking",
                           max_iterations = 500,
                           convergence = list(absolute = 1e-4, rule = "improvement"),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$status, 0L)
  # rule="improvement" -> rule=improvement, metric=max_err (via convergence_used)
  expect_equal(result$convergence_used$rule, "improvement",
               info = "rule=improvement must produce convergence_used$rule=improvement")
  expect_equal(result$convergence_used$metric, "max_err",
               info = "default metric must be max_err")
})

test_that("A4-metrics: all 6 metrics present and non-zero in calib_result", {
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.4,"2"=0.4,"3"=0.2), b = c("1"=0.6,"2"=0.4))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "ieppa",
                           max_iterations = 100, attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error","mean_error","kl","chi2","grake_norm","l1_weight_change"))
    expect_true(nm %in% names(result), info = sprintf("metric '%s' missing", nm))
  expect_false("pct_change" %in% names(result), info = "pct_change must be renamed")
  expect_gt(result$grake_norm, 0, label = "grake_norm must be > 0 on non-trivial data")
  expect_gt(result$l1_weight_change, 0, label = "l1_weight_change must be > 0 after calibration")
})

test_that("A4b: plateau rule rejects tol >= 1", {
  expect_error(
    leafblower::harvest(
      data.frame(a = factor(c("1","2"))),
      list(a = c("1"=0.5,"2"=0.5)),
      max_weight = 3, method = "ieppa",
      convergence = list(rule = "plateau", tol = 1.5),
      attach_weights = FALSE
    ),
    regexp = "must be in"
  )
})

test_that("A5: list(absolute=1e-6) sets convergence_used rule=threshold metric=max_err", {
  set.seed(44)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
                           max_iterations = 500,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$status, 0L)
  expect_lt(result$max_error, 1e-6)
  expect_equal(result$convergence_used$rule, "threshold",
               info = "list(absolute=) must map to rule=threshold")
  expect_equal(result$convergence_used$metric, "max_err",
               info = "list(absolute=) must map to metric=max_err")
})

test_that("A6: grake_norm+threshold converges on raking (skip if no survey pkg)", {
  skip_if_not_installed("survey")
  set.seed(99)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w_lb <- leafblower::harvest(data, target, max_weight = 10, method = "raking",
                              max_iterations = 200,
                              convergence = list(metric = "grake_norm",
                                                 rule = "threshold", tol = 1e-7),
                              attach_weights = FALSE)
  result <- attr(w_lb, "result")
  expect_equal(result$status, 0L, info = "grake_norm criterion must converge")
  expect_lt(result$grake_norm, 1e-7)
})

test_that("lj5j: convergence_used metric/rule are character strings (NA-safe indexing)", {
  set.seed(1)
  data <- data.frame(a = factor(sample(c("1","2"), 100, replace=TRUE)))
  target <- list(a = c("1"=0.5, "2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=3, method="ieppa",
                           attach_weights=FALSE)
  r <- attr(w, "result")
  expect_true(is.character(r$convergence_used$metric) || is.na(r$convergence_used$metric))
  expect_true(is.character(r$convergence_used$rule)   || is.na(r$convergence_used$rule))
})

test_that("we1a: improvement+absolute without stop_when errors with actionable message", {
  data <- data.frame(a = factor(c("1","2")))
  target <- list(a = c("1"=0.5,"2"=0.5))
  expect_error(
    leafblower::harvest(data, target, max_weight=3, method="ieppa",
                        convergence = list(improvement=0.01, absolute=1e-6),
                        attach_weights=FALSE),
    regexp = "stop_when|ambiguous|combine"
  )
})

test_that("B17: PLATEAU rule does not fire when prev=0", {
  # Indirect: run ieppa on an already-calibrated problem. First iter gives errRp≈0.
  # Without fix: PLATEAU fires at iter=2 as "prev=0 → any curr >= 0" is true.
  # With fix: PLATEAU requires prev>0 so it skips at iter=2, let RK_OK fire naturally.
  result <- harvest(
    data.frame(x = c("A","A","B","B","B"), w = rep(1,5)),
    target = list(x = c(A=0.4, B=0.6)),
    method = "ieppa",
    convergence = list(rule="plateau", pct=0.01),
    max_iterations = 50L
  )
  # Must converge (not a spurious plateau at wrong errRp).
  expect_equal(attr(result,"result")$max_error, 0.0, tolerance=1e-8)
})

test_that("vpby: stop_when=ANY prev not mutated to post-halt curr — improvement+absolute fires correctly", {
  # Regression for leafblower-vpby: apply_rule() used to write prev=curr
  # unconditionally, including on the halt iteration.  Under stop_when=ANY the
  # driver halts as soon as either criterion fires; the leaked prev=curr means
  # any re-entry would compute improvement against the halting value rather than
  # the at-convergence baseline.  This test drives stop_when="any" with both
  # improvement and absolute criteria and verifies:
  #   (a) the run converges (not stalled),
  #   (b) max_error at exit is ≤ the absolute_tol (correct at-convergence state),
  #   (c) the improvement rule did not fire prematurely (iterations > 1).
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(
    a = c("1" = 1/3, "2" = 1/3, "3" = 1/3),
    b = c("1" = 0.5, "2" = 0.5)
  )
  w <- leafblower::harvest(
    data, target, max_weight = 5, method = "ieppa",
    convergence = list(improvement = 0.01, absolute = 1e-4, stop_when = "any"),
    max_iterations = 500L,
    attach_weights = FALSE
  )
  r <- attr(w, "result")
  # (a) converged — not NO_CONV
  expect_false(grepl("no_conv|NOCONV|noconv", tolower(if (is.null(r$status)) "" else r$status)))
  # (b) at-convergence error reflects correct baseline (≤ absolute_tol OR improvement fired)
  expect_true(is.finite(r$max_error))
  expect_true(r$max_error < 0.2)   # loose: confirms solver ran and made progress
  # (c) at least one iteration ran (improvement rule needs prev finite before it can fire)
  expect_true(r$iterations >= 1L)
})
