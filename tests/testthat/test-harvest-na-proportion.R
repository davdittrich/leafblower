# CR-F1 (dtkn.1): add_na_proportion=TRUE must build a named NUMERIC target
# vector, not a list. The old c(lapply(...), list("NA"=)) produced a LIST, so
# compute_quality_metrics' margin_kl_one arithmetic hit "non-numeric argument to
# binary operator", was swallowed by the tryCatch at harvest.R:1113, and every
# add_na_proportion=TRUE run reported margin_kl = NA plus a spurious warning.
#
# NOTE: the *value* of margin_kl for an NA-injected margin is +Inf (not finite)
# because compute_quality_metrics excludes NA obs from the weighted marginal so
# the injected "NA" target bin has no mass — that is a SEPARATE bug tracked as
# dtkn.12 (CR-F1b). This test asserts only the dtkn.1 scope: the type error is
# gone (no "non-numeric" warning) and margin_kl is a real number, not NA.

test_that("CR-F1: add_na_proportion=TRUE produces numeric target, no non-numeric warning", {
  set.seed(7); n <- 400L
  df <- data.frame(
    a = factor(sample(c("x", "y", "z"), n, TRUE)),
    b = factor(sample(c("M", "F"), n, TRUE))
  )
  df$a[sample(n, 40L)] <- NA  # 10% NA in the target variable
  tgt <- list(a = c(x = 0.5, y = 0.3, z = 0.2), b = c(M = 0.6, F = 0.4))

  warns <- character(0)
  res <- withCallingHandlers(
    harvest(df, tgt, add_na_proportion = TRUE, attach_weights = FALSE),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  mkl <- attr(res, "result")$margin_kl

  # (a) the swallowed type error no longer fires
  expect_false(any(grepl("non-numeric|margin_kl:", warns)))
  # (b) margin_kl is a real number, not NA (dtkn.1 fixed the list->numeric type)
  expect_false(is.na(mkl))
})
