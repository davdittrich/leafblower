## test-calib-result-consolidation.R
## Regression tests for CalibResult base struct consolidation (ztid.4).
## Each solver must populate the shared CalibResult fields correctly.
library(leafblower)

make_data <- function(n = 50L) {
  set.seed(42L)
  wts <- runif(n, 0.8, 1.2)
  df  <- data.frame(x = factor(sample(c("a", "b"), n, TRUE)))
  pop <- list(x = c(a = 0.5, b = 0.5))
  list(df = df, pop = pop, wts = wts, n = n)
}

for (m in c("ieppa", "chebyshev", "raking", "sinkhorn", "greg",
            "greenkhorn", "logit")) {
  local({
    method <- m
    d <- make_data()
    test_that(paste("CalibResult fields populated for method:", method), {
      r <- harvest(d$df, d$pop,
                   method           = method,
                   max_iterations   = 300L,
                   base_weights     = d$wts)
      raw <- attr(r, "result")
      # iterations must be positive
      expect_true(raw$iterations > 0L,
                  label = paste(method, "iterations > 0"))
      # max_error must be finite
      expect_true(is.finite(raw$max_error),
                  label = paste(method, "max_error is finite"))
    })
  })
}
