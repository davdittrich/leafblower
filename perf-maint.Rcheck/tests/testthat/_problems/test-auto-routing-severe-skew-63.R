# Extracted from test-auto-routing-severe-skew.R:63

# setup ------------------------------------------------------------------------
library(testthat)
test_env <- simulate_test_env(package = "leafblower", path = "..")
attach(test_env, warn.conflicts = FALSE)

# prequel ----------------------------------------------------------------------
library(leafblower)
make_zero_compress_df <- function(seed = 1L, n = 300L, K = 5L, ncat = 5L) {
  set.seed(seed)
  cats <- letters[seq_len(ncat)]
  df <- as.data.frame(lapply(seq_len(K),
    function(k) factor(sample(cats, n, TRUE))),
    stringsAsFactors = FALSE)
  names(df) <- paste0("m", seq_len(K))
  df
}

# test -------------------------------------------------------------------------
df <- make_zero_compress_df(seed = 3L)
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
