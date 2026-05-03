library(leafblower)

# WH-g: AUTO routing target-skew gate (Epic-H).
# Rule: K>=5 + M_cell/n>0.9 + target_skew>5 -> ieppa+sraa (else newton_kl).
# Fixtures use K=5 with 5-cat margins and n=300 to force M_cell/n ~ 1.0
# (5^5 = 3125 possible cells, ~300 unique tuples in random sample).

make_zero_compress_df <- function(seed = 1L, n = 300L, K = 5L, ncat = 5L) {
  set.seed(seed)
  cats <- letters[seq_len(ncat)]
  df <- as.data.frame(lapply(seq_len(K),
    function(k) factor(sample(cats, n, TRUE))),
    stringsAsFactors = FALSE)
  names(df) <- paste0("m", seq_len(K))
  df
}

test_that("WH-g: AUTO routes severe-skew K>=5 to ieppa", {
  df <- make_zero_compress_df(seed = 1L)
  # Severe skew per margin: max/min = 0.7/0.075 = 9.33 > 5
  tgt <- lapply(df, function(f) {
    p <- c(0.7, 0.075, 0.075, 0.075, 0.075)
    setNames(p, levels(f))
  })
  r <- suppressWarnings(harvest(df, tgt, method = "auto",
    max_weight = 10, min_weight = 0,
    max_iterations = 200L, attach_weights = FALSE, verbose = 0L))
  expect_equal(attr(r, "algorithm"), "ieppa",
    label = sprintf("WH-g severe-skew K=5: expected ieppa got %s",
                    attr(r, "algorithm")))
})

test_that("WH-g: AUTO keeps moderate-skew K>=5 on newton_kl", {
  df <- make_zero_compress_df(seed = 2L)
  # Moderate skew per margin: max/min = 0.4/0.1 = 4 <= 5
  tgt <- lapply(df, function(f) {
    p <- c(0.4, 0.2, 0.2, 0.1, 0.1)
    setNames(p, levels(f))
  })
  r <- suppressWarnings(harvest(df, tgt, method = "auto",
    max_weight = 10, min_weight = 0,
    max_iterations = 100L, attach_weights = FALSE, verbose = 0L))
  expect_equal(attr(r, "algorithm"), "newton_kl",
    label = sprintf("WH-g moderate-skew K=5: expected newton_kl got %s",
                    attr(r, "algorithm")))
})

test_that("WH-g: AUTO with min_target=0 takes severe-skew branch (no div-by-zero)", {
  df <- make_zero_compress_df(seed = 3L)
  # min target == 0: floor at 1e-12 forces target_skew large -> severe-skew
  # branch. ieppa may NOCONV on infeasible target; auto-fallback -> lbfgsb.
  # Either path is acceptable; the test checks NO CRASH and NOT-newton_kl.
  tgt <- lapply(df, function(f) {
    p <- c(0.55, 0.30, 0.10, 0.05, 0.0)
    setNames(p, levels(f))
  })
  r <- suppressWarnings(harvest(df, tgt, method = "auto",
    max_weight = 10, min_weight = 0,
    max_iterations = 300L, attach_weights = FALSE, verbose = 0L))
  alg <- attr(r, "algorithm")
  expect_true(alg %in% c("ieppa", "lbfgsb"),
    label = sprintf("WH-g zero-min-target: must not crash and not route newton_kl; got %s",
                    alg))
})
