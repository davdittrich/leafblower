# CR-D4 (leafblower-j7x8.4): raking and oris must run the shared capacity/
# negativity pre-entry feasibility check that the 6 full-setup solvers already run.
# A capacity-infeasible bound set (sum(L_c) > n) must return an INFEAS error up
# front (0 iters) instead of burning the full iteration budget and returning
# garbage weights with only a stall warning.

test_that("raking rejects capacity-infeasible bounds up front (CR-D4)", {
  set.seed(1); n <- 300
  df <- data.frame(a = factor(sample(c("x", "y"), n, TRUE)),
                   b = factor(sample(c("p", "q"), n, TRUE)),
                   c = factor(sample(c("m", "o"), n, TRUE)))
  tgt <- list(a = c(x = 0.5, y = 0.5), b = c(p = 0.5, q = 0.5), c = c(m = 0.5, o = 0.5))
  # min_weight=2 with mean-1 base ⇒ sum(L_c) = 2n > n ⇒ structurally infeasible.
  expect_error(
    suppressWarnings(harvest(df, tgt, method = "raking",
                             min_weight = 2, max_weight = 5, attach_weights = FALSE)),
    regexp = "iters|infeasible|bound|capacity|error"
  )
})

test_that("oris rejects capacity-infeasible bounds up front (CR-D4)", {
  set.seed(1); n <- 300
  df <- data.frame(a = factor(sample(c("x", "y"), n, TRUE)),
                   b = factor(sample(c("p", "q"), n, TRUE)),
                   c = factor(sample(c("m", "o"), n, TRUE)))
  tgt <- list(a = c(x = 0.5, y = 0.5), b = c(p = 0.5, q = 0.5), c = c(m = 0.5, o = 0.5))
  expect_error(
    suppressWarnings(harvest(df, tgt, method = "oris",
                             min_weight = 2, max_weight = 5, attach_weights = FALSE)),
    regexp = "iters|infeasible|bound|capacity|error"
  )
})

test_that("raking/oris happy path still converges after the pre-entry check (CR-D4)", {
  set.seed(2); n <- 300
  df <- data.frame(a = factor(sample(c("x", "y"), n, TRUE)),
                   b = factor(sample(c("p", "q"), n, TRUE)))
  tgt <- list(a = c(x = 0.5, y = 0.5), b = c(p = 0.5, q = 0.5))
  rk <- attr(suppressWarnings(harvest(df, tgt, method = "raking", max_iter = 200)), "result")
  or <- attr(suppressWarnings(harvest(df, tgt, method = "oris",   max_iter = 200)), "result")
  expect_equal(rk$status, 0L)   # RK_OK — feasible problem unaffected by the guard
  expect_equal(or$status, 0L)
})
