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
  # l1_weight_change must be computed (non-negative) and within tolerance
  # (may be exactly 0 when solver converges before the first pct check interval)
  expect_gte(result$l1_weight_change, 0)
  expect_lt(result$l1_weight_change, 0.001 * 1.5)
  expect_lt(result$max_error, 1e-3)
  expect_equal(result$status, 0L)
  # Verify pct_change IS computed: use a barely-converged run (max_iterations=10)
  # where the solver stops before pct threshold is met; l1_weight_change must be > 0
  w10 <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    max_iterations = 10, attach_weights = FALSE)
  r10 <- attr(w10, "result")
  expect_gt(r10$l1_weight_change, 0)
})

test_that("A8a: criterion='mean_err' actively stops solver", {
  set.seed(43)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    convergence = list(absolute = 1e-4, criterion = "mean_err"),
    max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  # mean_error must be computed (non-zero) and within threshold
  expect_gt(result$mean_error, 0)
  expect_lt(result$mean_error, 1e-4)
  expect_equal(result$status, 0L)
})

test_that("A8b: criterion='kl' actively stops solver", {
  set.seed(44)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    convergence = list(absolute = 1e-6, criterion = "kl"),
    max_iterations = 500, attach_weights = FALSE)
  result <- attr(w, "result")
  # kl must be computed (non-zero) and within threshold
  expect_gt(result$kl, 0)
  expect_lt(result$kl, 1e-6)
  expect_equal(result$status, 0L)
})

test_that("A8c: criterion='chi2' actively stops solver", {
  set.seed(45)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  # chi2 scales with n; threshold n-scaled (~1e-3 * 2000 = 2)
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    convergence = list(absolute = 2.0, criterion = "chi2"),
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

test_that("hawe: iEPPA warns on PCT stall with large max_error", {
  # Contradictory targets: var1 wants A=95% but var2 structure makes that impossible
  n <- 400
  var1 <- factor(rep(c("A","B"), each = n/2))
  var2 <- factor(rep(c("2","1"), each = n/2))  # A maps to "2", B maps to "1"
  data <- data.frame(var1 = var1, var2 = var2)
  target <- list(
    var1 = c(A = 0.95, B = 0.05),
    var2 = c("1" = 0.95, "2" = 0.05)
  )
  expect_warning(
    leafblower::harvest(data, target, max_weight = 1.5, method = "ieppa",
                        max_iterations = 300,
                        convergence = list(pct = 0.001),
                        attach_weights = FALSE),
    regexp = "PCT convergence stall"
  )
})

test_that("A1: default convergence (max_err+improvement) converges smooth synthetic", {
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
  # default: rule=IMPROVEMENT (1), metric=max_err (0)
  expect_equal(result$convergence_rule, 1L,
               info = "default rule must be IMPROVEMENT (1)")
  expect_equal(result$convergence_metric, 0L,
               info = "default metric must be max_err (0)")
})

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
  if (result$status == 1L) {
    expect_lt(result$best_error, 0.9 * result$max_error)
  }
})

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
  expect_lt(result$l1_weight_change, 0.001)
  # pct key -> metric=l1_weight (5), rule=plateau (2) — autumn/anesrake compatible
  expect_equal(result$convergence_metric, 5L,
               info = "pct key must select l1_weight metric (5)")
  expect_equal(result$convergence_rule, 2L,
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
  # rule="improvement" -> rule=improvement (1), metric=max_err (0)
  expect_equal(result$convergence_rule, 1L,
               info = "rule=improvement must produce convergence_rule=1")
  expect_equal(result$convergence_metric, 0L,
               info = "default metric must be max_err (0)")
})
