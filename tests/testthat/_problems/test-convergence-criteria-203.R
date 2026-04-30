# Extracted from test-convergence-criteria.R:203

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "leafblower", path = "..")
attach(test_env, warn.conflicts = FALSE)

# test -------------------------------------------------------------------------
n <- 400
var1 <- factor(rep(c("A","B"), each = n/2))
var2 <- factor(rep(c("2","1"), each = n/2))
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
