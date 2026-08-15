# Tests for RVAL.1, RVAL.2, RVAL.3, META.2

# ---------------------------------------------------------------------------
# RVAL.1: validate target sums to 1 before add_na_proportion rescale
# ---------------------------------------------------------------------------
test_that("RVAL.1: add_na_proportion errors when target does not sum to 1", {
  set.seed(1)
  n <- 100L
  df <- data.frame(x = c(rep("a", 50), rep("b", 40), rep(NA, 10)))
  # Partial target: sums to 0.7, not 1
  tgt <- list(x = c(a = 0.4, b = 0.3))
  expect_error(
    harvest(df, tgt, add_na_proportion = TRUE,
            convergence = list(absolute = 1e-4)),
    regexp = "target.*x.*sum|sum.*1.*x",
    ignore.case = TRUE
  )
})

test_that("RVAL.1: add_na_proportion accepts valid target summing to 1", {
  set.seed(1)
  n <- 100L
  df <- data.frame(x = c(rep("a", 50), rep("b", 40), rep(NA, 10)))
  # Valid target: sums to 1
  tgt <- list(x = c(a = 0.55, b = 0.45))
  expect_no_error(
    harvest(df, tgt, add_na_proportion = TRUE,
            convergence = list(absolute = 1e-4))
  )
})

# ---------------------------------------------------------------------------
# RVAL.2: warn on unknown harvest() ... args
# ---------------------------------------------------------------------------
test_that("RVAL.2: unknown ... arg emits warning listing the arg name", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a", "b"), 100, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  expect_warning(
    harvest(df, tgt, convergence = list(absolute = 1e-4),
            typo_arg = 99),
    regexp = "typo_arg",
    ignore.case = TRUE
  )
})

test_that("RVAL.2: normal call with no ... does not emit unknown-arg warning", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a", "b"), 100, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  # Suppress unrelated warnings (budget, sparse, etc.)
  warns <- character(0)
  withCallingHandlers(
    harvest(df, tgt, convergence = list(absolute = 1e-4)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  unknown_warns <- warns[grepl("unknown.*arg|unrecognized|ignored", warns, ignore.case = TRUE)]
  expect_length(unknown_warns, 0)
})

# ---------------------------------------------------------------------------
# RVAL.3: warn on OOV observations in group_ids_r
# ---------------------------------------------------------------------------
test_that("RVAL.3: unlisted level in data emits OOV warning", {
  set.seed(1)
  # "c" is in data but not in target
  df  <- data.frame(x = c("a", "b", "c", "a", "b"))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  # Under testthat 3e, expect_warning() only muffles the ONE warning matching
  # its regexp and lets harvest()'s incidental sparse-category and
  # fixed-point-convergence diagnostics propagate -- edition 2 swallowed all
  # of them (D-11, phase 01-04). Muffle ONLY those two named incidental
  # patterns so edition 3's stricter unmatched-warning signal still catches
  # an unrelated regression in this call.
  withCallingHandlers(
    expect_warning(
      harvest(df, tgt, convergence = list(absolute = 1e-4)),
      regexp = "out.of.vocabulary|OOV|unlisted|unknown level",
      ignore.case = TRUE
    ),
    warning = function(w) {
      if (grepl("sparse categories detected|budget exhausted|did not converge|plateau|fixed point at",
                 conditionMessage(w)))
        invokeRestart("muffleWarning")
    }
  )
})

test_that("RVAL.3: NA obs on an NA-margin does NOT trigger OOV warning", {
  set.seed(1)
  df <- data.frame(x = c(rep("a", 50), rep("b", 40), rep(NA, 10)))
  tgt <- list(x = c(a = 0.55, b = 0.45))
  warns <- character(0)
  withCallingHandlers(
    harvest(df, tgt, add_na_proportion = TRUE,
            convergence = list(absolute = 1e-4)),
    warning = function(w) {
      warns <<- c(warns, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  oov_warns <- warns[grepl("out.of.vocabulary|OOV|unlisted|unknown level", warns,
                           ignore.case = TRUE)]
  expect_length(oov_warns, 0)
})

# ---------------------------------------------------------------------------
# META.2: pct alias in metric_int — resolve earlier, remove duplicate
# ---------------------------------------------------------------------------
test_that("META.2: convergence metric='pct' selects l1_weight (integer 5L)", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a", "b"), 200, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  # Should complete without error; pct maps to l1_weight=5L
  expect_no_error(
    harvest(df, tgt,
            convergence = list(metric = "pct", absolute = 1e-4))
  )
})

test_that("META.2: convergence list(pct=0.005) still works (pct shorthand)", {
  set.seed(1)
  df  <- data.frame(x = factor(sample(c("a", "b"), 200, replace = TRUE)))
  tgt <- list(x = c(a = 0.5, b = 0.5))
  expect_no_error(
    harvest(df, tgt,
            convergence = list(pct = 0.005))
  )
})

test_that("dtkn.6: plain NAs do NOT trigger the OOV warning (add_na_proportion=FALSE) (CR-F6)", {
  set.seed(1L); n <- 200L
  a <- factor(sample(c("x", "y", "z"), n, TRUE)); a[sample(n, 20L)] <- NA
  df <- data.frame(a = a); tg <- list(a = c(x = 0.34, y = 0.33, z = 0.33))
  # 20 plain NAs (normal missingness) must NOT warn about out-of-vocabulary.
  expect_no_warning(suppressMessages(harvest(df, tg, attach_weights = FALSE)))
})

test_that("dtkn.6: GENUINE out-of-vocabulary values still warn (CR-F6)", {
  # A factor level present in data but absent from the target IS a vocabulary problem.
  df <- data.frame(a = factor(c(rep("x", 90), rep("y", 90), rep("zzz", 20))))
  tg <- list(a = c(x = 0.5, y = 0.5))  # 'zzz' not in target
  expect_warning(suppressMessages(harvest(df, tg, attach_weights = FALSE)),
                 "out-of-vocabulary")
})
