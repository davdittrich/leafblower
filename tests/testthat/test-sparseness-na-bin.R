## Tests for compute_sparseness_diag NA-bin counting — 4ihf.6
## Regression guard for the true-NA-blind bug: tabulate(match(col, lvs)) dropped
## true-NA observations (match(NA, lvs) -> NA_integer_), so n_kj("NA") was 0 on
## EVERY add_na_proportion run, mis-flagging the NA bin sparse AND diverging from
## the solver (which conflates true-NA + literal-"NA" into one "NA" group) and
## from Python (pd.isna().sum()). The fix mirrors the solver encoding for
## injected margins: as.character(col); key[is.na(key)] <- "NA".

library(leafblower)

# Fixture with BOTH true-NA AND a literal-string-"NA" category in one margin.
.make_na_fixture <- function() {
  # 3 true-NA + 2 literal-"NA" + 5 "x" + 1 "y"  (n = 11)
  g <- c("x", "x", "x", "x", "x", "y", "NA", "NA", NA, NA, NA)
  list(
    df     = data.frame(g = g, stringsAsFactors = FALSE),
    # target before injection sums to 1; the solver path injects the "NA" bin.
    target = list(g = c(x = 0.7, y = 0.3))
  )
}

test_that("compute_sparseness_diag: conflates true-NA + literal-'NA' for injected margin", {
  fx  <- .make_na_fixture()
  df  <- fx$df
  csd <- leafblower:::compute_sparseness_diag

  # Build the post-injection target exactly as harvest() does (yaye path).
  na_frac <- mean(is.na(df$g))                       # 3/11
  tgt_inj <- list(g = c(
    lapply(fx$target$g, function(t) t * (1 - na_frac)),
    list("NA" = na_frac)
  ))
  tgt_inj$g <- unlist(tgt_inj$g)

  conflated <- sum(is.na(df$g)) + sum(!is.na(df$g) & as.character(df$g) == "NA")
  expect_equal(conflated, 5L)                        # 3 true-NA + 2 literal-"NA"

  # Low obs_threshold so the NA entry is reported with its count regardless of sparsity.
  diag_inj <- csd(tgt_inj, df, cat_threshold = 0.99, obs_threshold = 100L,
                  na_margins = "g")
  na_entry <- Filter(function(e) e$level == "NA", diag_inj$g)
  expect_length(na_entry, 1L)
  na_count <- na_entry[[1]]$n_kj

  # The fix: conflated count, NOT the old true-NA-blind 0.
  expect_equal(na_count, conflated)                  # 5
  expect_false(na_count == 0L)                       # bug was n_kj("NA") == 0

  # Matches the solver's NA-group size: group_ids_r else-branch (harvest.R:475-477)
  # maps both true-NA and literal-"NA" to the "NA" target index.
  key <- as.character(df$g)
  key[is.na(key)] <- "NA"
  solver_na_group <- sum(match(key, names(tgt_inj$g)) == which(names(tgt_inj$g) == "NA"))
  expect_equal(na_count, solver_na_group)            # 5
})

test_that("compute_sparseness_diag: NA bin not falsely sparse when above obs_threshold", {
  # 40 true-NA + 5 literal-"NA" = 45 > obs_threshold(30); 100 "x".
  g <- c(rep("x", 100), rep("NA", 5), rep(NA, 40))
  df <- data.frame(g = g, stringsAsFactors = FALSE)
  na_frac <- mean(is.na(df$g))
  tgt <- list(g = c(x = (1 - na_frac), "NA" = na_frac))
  csd <- leafblower:::compute_sparseness_diag

  diag <- csd(tgt, df, cat_threshold = 0.01, obs_threshold = 30L, na_margins = "g")
  sparse_levels <- vapply(diag$g, function(e) e$level, character(1))
  expect_false("NA" %in% sparse_levels)              # 45 obs — not sparse
})

test_that("compute_sparseness_diag: non-injected margins are byte-unchanged (no na_margins)", {
  # Margin with true-NA but NOT injected: keeps the original tabulate(match(col, lvs))
  # behaviour — true-NA is NOT counted under any level (match(NA, lvs) -> NA).
  g <- c(rep("x", 5), rep("y", 5), rep(NA, 5))
  df <- data.frame(g = g, stringsAsFactors = FALSE)
  tgt <- list(g = c(x = 0.5, y = 0.5))
  csd <- leafblower:::compute_sparseness_diag

  # High obs_threshold so every level is reported with its count.
  diag_default <- csd(tgt, df, cat_threshold = 0.01, obs_threshold = 100L)
  diag_empty_na <- csd(tgt, df, cat_threshold = 0.01, obs_threshold = 100L,
                       na_margins = character(0))
  expect_identical(diag_default, diag_empty_na)      # default == explicit empty
  # x and y each have 5 obs; no "NA" level in target -> true-NA simply uncounted.
  x_entry <- Filter(function(e) e$level == "x", diag_default$g)
  expect_equal(x_entry[[1]]$n_kj, 5L)
})
