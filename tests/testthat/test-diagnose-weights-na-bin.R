## NABIN.2 — diagnose_weights must count the explicit "NA" bin (leafblower-4ihf.2)
## Same root cause as get_current_miss: `not_na & (col_char == "NA")` is always
## FALSE, so prop_original / prop_weighted for the injected "NA" bin reported 0.

library(leafblower)

test_that("diagnose_weights: 'NA' bin prop_original ~ na_frac, not 0", {
  set.seed(13L)
  n <- 500L
  na_frac <- 0.20
  na_count <- round(n * na_frac)
  g <- sample(c("a", "b", "c"), n, replace = TRUE)
  na_idx <- sample(n, na_count)
  g[na_idx] <- NA
  df <- data.frame(g = g, stringsAsFactors = FALSE)

  # Target mirrors harvest(add_na_proportion = TRUE).
  target <- list(g = c(
    a    = 0.40 * (1 - na_frac),
    b    = 0.35 * (1 - na_frac),
    c    = 0.25 * (1 - na_frac),
    `NA` = na_frac
  ))

  # Uniform weights => prop_weighted == prop_original for every level.
  w <- rep(1.0, n)
  d <- diagnose_weights(df, target, w)

  na_row <- d[d$variable == "g" & d$level == "NA", ]
  expect_equal(nrow(na_row), 1L)

  # prop_original of the NA bin == na_count/n (all-obs denominator), not 0.
  expect_equal(na_row$prop_original, na_count / n, tolerance = 1e-12)
  expect_equal(na_row$prop_weighted, na_count / n, tolerance = 1e-12)

  # Shares sum to 1 across all bins (all-obs denominator).
  expect_equal(sum(d$prop_original), 1.0, tolerance = 1e-12)
  expect_equal(sum(d$prop_weighted), 1.0, tolerance = 1e-12)
})

test_that("diagnose_weights: NO 'NA' bin but data HAS NA -> shares over non-NA, error_weighted ~ 0", {
  ## Regression guard (NABIN): when the target has no explicit "NA" level but
  ## the data contains NA, the named-level denominator must EXCLUDE NA obs.
  ## With an all-obs denominator the named shares sum to (1 - na_frac), so
  ## error_weighted would be ~ -na_frac*target on every level even for
  ## perfectly calibrated weights. The fix restores the non-NA denominator.
  set.seed(29L)
  n <- 500L
  na_frac <- 0.20
  na_count <- round(n * na_frac)
  g <- sample(c("a", "b", "c"), n, replace = TRUE)
  na_idx <- sample(n, na_count)
  g[na_idx] <- NA
  df <- data.frame(g = g, stringsAsFactors = FALSE)

  # Target over the NON-NA observed shares (no "NA" bin). Set the target equal
  # to the observed non-NA shares so well-calibrated (here: uniform) weights
  # give error_weighted == 0.
  n_nonna <- sum(!is.na(g))
  obs_a <- sum(g == "a", na.rm = TRUE) / n_nonna
  obs_b <- sum(g == "b", na.rm = TRUE) / n_nonna
  obs_c <- sum(g == "c", na.rm = TRUE) / n_nonna
  target <- list(g = c(a = obs_a, b = obs_b, c = obs_c))

  w <- rep(1.0, n)
  d <- diagnose_weights(df, target, w)

  # Named-level shares sum to 1 over NON-NA obs (NOT 1 - na_frac).
  expect_equal(sum(d$prop_original), 1.0, tolerance = 1e-12)
  expect_equal(sum(d$prop_weighted), 1.0, tolerance = 1e-12)

  # error_weighted ~ 0 on EVERY level (the regression produced ~ -na_frac*tgt).
  expect_true(all(abs(d$error_weighted) < 1e-12))
  expect_true(all(abs(d$error_original) < 1e-12))
})
