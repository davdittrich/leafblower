# Extracted from test-calibration-solvers.R:253

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "leafblower", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(5)
n <- 400
data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
target <- list(a=c("1"=0.4,"2"=0.4,"3"=0.2), b=c("1"=0.6,"2"=0.4))
w_greg <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                                method="greg", attach_weights=FALSE)
w_rake <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                                method="raking", max_iterations=500, attach_weights=FALSE)
r_greg <- attr(w_greg, "result")
r_rake <- attr(w_rake, "result")
expect_equal(r_greg$status, 0L, info="greg must converge")
expect_equal(r_rake$status, 0L, info="raking must converge")
expect_lte(r_greg$chi2, r_rake$chi2 + 1e-6,
             label="greg chi2 <= raking chi2")
