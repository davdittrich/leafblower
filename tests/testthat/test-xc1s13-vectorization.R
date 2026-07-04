# ──────────────────────────────────────────────────────────────────────────────
# xc1s.13(f)/(g): permanent invariants for the single-pass rewrites of
# weighted_pct (rowsum), diagnose_weights, and get_current_miss. Guards the
# NA/empty-level/NA-bin edge cases that byte-identity depends on.
# ──────────────────────────────────────────────────────────────────────────────

test_that("weighted_pct: empty level→0, NA excluded, in-vocab sums to 1 (xc1s.13f)", {
  x <- factor(c("a", "a", "b", NA, "a"), levels = c("a", "b", "z"))  # 'z' empty, one NA
  w <- c(1, 2, 3, 4, 5)
  p <- weighted_pct(x, w)
  expect_named(p, c("a", "b", "z"))
  expect_equal(unname(p[["z"]]), 0)          # empty level → 0
  # non-NA weight total = 1+2+3+5 = 11; a = 1+2+5 = 8; b = 3
  expect_equal(unname(p[["a"]]), 8 / 11)
  expect_equal(unname(p[["b"]]), 3 / 11)
  expect_equal(sum(p), 1)                    # NA row excluded from denominator
})

test_that("weighted_pct: zero total weight → all-zero (xc1s.13f)", {
  x <- factor(c("a", "b"), levels = c("a", "b"))
  expect_equal(unname(weighted_pct(x, c(0, 0))), c(0, 0))
})

test_that("diagnose_weights: NA bin conflates true-missing + literal 'NA' (xc1s.13g)", {
  df <- data.frame(a = c("x", "x", NA, "NA", "y"), stringsAsFactors = FALSE)  # 1 true-NA, 1 literal
  w  <- c(1, 1, 1, 1, 1)
  tg <- list(a = c(x = .4, y = .2, "NA" = .4))       # explicit NA bin
  d  <- diagnose_weights(df, tg, w)
  get <- function(lv) d[d$variable == "a" & d$level == lv, ]
  expect_equal(get("NA")$prop_original, 2 / 5)       # true-missing + literal-"NA" conflated
  expect_equal(get("x")$prop_original,  2 / 5)       # literal-"NA" NOT counted as its own level
  expect_equal(get("y")$prop_original,  1 / 5)
  expect_equal(get("NA")$prop_weighted, 2 / 5)       # equal weights
})

test_that("diagnose_weights: empty target level → 0 proportion (xc1s.13g)", {
  df <- data.frame(b = factor(c("m", "n", "m"), levels = c("m", "n", "q")), stringsAsFactors = FALSE)
  tg <- list(b = c(m = .5, n = .3, q = .2))          # 'q' never observed
  d  <- diagnose_weights(df, tg, runif(3))
  qr <- d[d$variable == "b" & d$level == "q", ]
  expect_equal(qr$prop_original, 0)
  expect_equal(qr$prop_weighted, 0)
})

test_that("get_current_miss: weighted max miss incl. NA bin (xc1s.13g)", {
  df <- data.frame(a = c("x", "x", NA, "NA", "y"), stringsAsFactors = FALSE)
  w  <- c(1, 1, 1, 1, 1)
  tg <- list(a = c(x = .4, y = .2, "NA" = .4))
  # get_current_miss returns the per-margin max miss, named by margin.
  # weighted props: x=2/5=.4 (miss 0), y=1/5=.2 (miss 0), NA=2/5=.4 (miss 0) → max 0
  expect_equal(unname(get_current_miss(df, tg, w)), 0)
  # skew a target to force a known miss
  tg2 <- list(a = c(x = .9, y = .05, "NA" = .05))
  expect_equal(unname(get_current_miss(df, tg2, w)), max(abs(c(.4 - .9, .2 - .05, .4 - .05))))
})
