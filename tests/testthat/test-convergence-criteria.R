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
  # pct_change must be computed (non-negative) and within tolerance
  # (may be exactly 0 when solver converges before the first pct check interval)
  expect_gte(result$pct_change, 0)
  expect_lt(result$pct_change, 0.001 * 1.5)
  expect_lt(result$max_error, 1e-3)
  expect_equal(result$status, 0L)
  # Verify pct_change IS computed: use a barely-converged run (max_iterations=10)
  # where the solver stops before pct threshold is met; pct_change must be > 0
  w10 <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
    max_iterations = 10, attach_weights = FALSE)
  r10 <- attr(w10, "result")
  expect_gt(r10$pct_change, 0)
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
