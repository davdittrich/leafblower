## eb79.9 — contract test for .encode_na_bin_mask(), the shared internal
## helper factored out of diagnose_weights.R and current_miss.R's identical
## "conflate true-NA with literal string 'NA'" mask blocks.
## Internal (non-exported) helper: bind via ::: so the bare calls below resolve
## under the installed-package namespace (devtools::test load_all also works).
.encode_na_bin_mask <- leafblower:::.encode_na_bin_mask

test_that(".encode_na_bin_mask: target_has_na = TRUE conflates true-NA and literal 'NA'", {
  col <- c("a", NA, "NA", "b", NA, "b")
  mask <- .encode_na_bin_mask(col, TRUE)
  # rows 2, 3, 5 are true-NA or literal "NA"; rows 1, 4, 6 are normal levels.
  expect_equal(mask, c(FALSE, TRUE, TRUE, FALSE, TRUE, FALSE))
})

test_that(".encode_na_bin_mask: target_has_na = FALSE returns plain is.na(col)", {
  col <- c("a", NA, "NA", "b", NA, "b")
  mask <- .encode_na_bin_mask(col, FALSE)
  # No conflation: only genuinely-missing rows are TRUE; literal "NA" string
  # (row 3) is treated as a normal (non-missing) level value.
  expect_equal(mask, is.na(col))
  expect_equal(mask, c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE))
})

test_that(".encode_na_bin_mask: all-normal column has no NA-bin hits either way", {
  col <- c("a", "b", "c")
  expect_equal(.encode_na_bin_mask(col, TRUE), c(FALSE, FALSE, FALSE))
  expect_equal(.encode_na_bin_mask(col, FALSE), c(FALSE, FALSE, FALSE))
})

test_that(".encode_na_bin_mask: numeric column with NA (non-character input)", {
  col <- c(1, NA, 3, NA)
  # as.character(NA_real_) is NA_character_, never the literal "NA" string,
  # so numeric columns never trigger the literal-string branch either way.
  expect_equal(.encode_na_bin_mask(col, TRUE), c(FALSE, TRUE, FALSE, TRUE))
  expect_equal(.encode_na_bin_mask(col, FALSE), c(FALSE, TRUE, FALSE, TRUE))
})
