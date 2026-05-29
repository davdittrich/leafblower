## NABIN.1 — get_current_miss must count the explicit "NA" bin (leafblower-4ihf.1)
## as.character(NA) is NA_character_, and `!is.na(col) & as.character(col)=="NA"`
## is ALWAYS FALSE, so the injected "NA" target bin previously reported 0 share.

library(leafblower)

test_that("get_current_miss: 'NA' bin observed share is counted (~ na_frac), not 0", {
  set.seed(11L)
  n <- 500L
  na_frac <- 0.20
  na_count <- round(n * na_frac)
  g <- sample(c("a", "b", "c"), n, replace = TRUE)
  na_idx <- sample(n, na_count)
  g[na_idx] <- NA
  df <- data.frame(g = g, stringsAsFactors = FALSE)

  # Construct a target whose NON-NA levels exactly match the observed all-obs
  # shares, so the ONLY source of miss is the "NA" bin. With the bug the "NA"
  # bin observed share is 0 (miss == na_frac); with the fix it is na_frac so
  # we can set the NA target to that and obtain miss ~ 0.
  obs_a <- sum(g == "a", na.rm = TRUE) / n
  obs_b <- sum(g == "b", na.rm = TRUE) / n
  obs_c <- sum(g == "c", na.rm = TRUE) / n
  obs_na <- na_count / n

  # Bug-exposing target: set NA target to its true observed share. If the NA
  # bin is ignored, get_current_miss reports |0 - obs_na| = obs_na (~0.20).
  target <- list(g = c(a = obs_a, b = obs_b, c = obs_c, `NA` = obs_na))

  w <- rep(1.0, n)
  miss <- get_current_miss(df, target, w)

  # With the fix every bin (including "NA") matches its observed share exactly.
  expect_lt(unname(miss[["g"]]), 1e-12)

  # Independent guard: the NA-bin observed share must equal na_frac, not 0.
  expect_equal(obs_na, na_frac, tolerance = 1e-12)
})

test_that("get_current_miss: NO 'NA' bin but data HAS NA -> miss over non-NA, ~ 0", {
  ## Regression guard (NABIN): with no explicit "NA" target level the
  ## denominator must EXCLUDE NA obs. An all-obs denominator deflates the
  ## named-level shares by (1 - na_frac), so a perfectly calibrated target
  ## would report miss ~ na_frac*target on every level. The fix restores the
  ## non-NA denominator => miss ~ 0.
  set.seed(31L)
  n <- 500L
  na_frac <- 0.20
  na_count <- round(n * na_frac)
  g <- sample(c("a", "b", "c"), n, replace = TRUE)
  na_idx <- sample(n, na_count)
  g[na_idx] <- NA
  df <- data.frame(g = g, stringsAsFactors = FALSE)

  # Target = observed NON-NA shares (no "NA" bin).
  n_nonna <- sum(!is.na(g))
  obs_a <- sum(g == "a", na.rm = TRUE) / n_nonna
  obs_b <- sum(g == "b", na.rm = TRUE) / n_nonna
  obs_c <- sum(g == "c", na.rm = TRUE) / n_nonna
  target <- list(g = c(a = obs_a, b = obs_b, c = obs_c))

  w <- rep(1.0, n)
  miss <- get_current_miss(df, target, w)

  # With the non-NA denominator the observed shares match the target exactly.
  expect_lt(unname(miss[["g"]]), 1e-12)
})
