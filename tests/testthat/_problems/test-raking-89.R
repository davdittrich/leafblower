# Extracted from test-raking.R:89

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "leafblower", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
set.seed(91)
n <- 1000
df <- data.frame(cat = sample(c("A", "B"), n, replace = TRUE, prob = c(0.05, 0.95)))
tgt <- list(cat = c(A = 0.95, B = 0.05))
t0 <- Sys.time()
msgs <- capture.output(
    res <- suppressWarnings(harvest(df, tgt, method = "raking",
                                     max_weight = 1.2,
                                     max_iterations = 500,
                                     verbose = 1L,
                                     convergence = list(absolute = 1e-6))),
    type = "output"
  )
elapsed <- as.numeric(Sys.time() - t0, units = "secs")
expect_lt(elapsed, 5)
probe <- paste(msgs, collapse = "\n")
expect_match(probe, "errRp stalled for [0-9]+ consecutive checks",
               info = paste("expected descent-monitor message; got:", probe))
