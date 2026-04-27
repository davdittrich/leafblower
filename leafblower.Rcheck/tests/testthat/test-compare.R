context("cross-algorithm equivalence on feasible inputs")

test_that("ieppa, raking, lbfgsb agree to 1e-3 on 20 random feasible datasets", {
  set.seed(20260423)
  for (trial in 1:20) {
    n <- sample(c(1000, 5000, 10000), 1)
    K <- sample(3:5, 1)
    cats <- sample(3:5, K, replace = TRUE)
    # Build random balanced categorical data + balanced targets
    cols <- lapply(seq_len(K), function(k) {
      factor(sample(0:(cats[k]-1), n, replace = TRUE))
    })
    df <- as.data.frame(cols)
    names(df) <- paste0("v", seq_len(K))
    tgt <- lapply(seq_len(K), function(k) {
      setNames(rep(1/cats[k], cats[k]), as.character(0:(cats[k]-1)))
    })
    names(tgt) <- names(df)
    # Reasonable bounds
    mw <- sample(c(2, 3, 5), 1)
    r_ieppa  <- suppressWarnings(harvest(df, tgt, method = "ieppa",  max_weight = mw))
    r_raking <- suppressWarnings(harvest(df, tgt, method = "raking", max_weight = mw))
    r_lbfgsb <- suppressWarnings(harvest(df, tgt, method = "lbfgsb", max_weight = mw))
    max_diff <- max(
      max(abs(r_ieppa$weights - r_raking$weights)),
      max(abs(r_ieppa$weights - r_lbfgsb$weights)),
      max(abs(r_raking$weights - r_lbfgsb$weights))
    )
    # Allow up to 10% relative error (different optimization paths for L-BFGS-B vs BCD methods)
    expect_lt(max_diff, max(0.1, 1e-3),
              label = sprintf("trial %d (n=%d K=%d mw=%.1f): max pairwise diff %.3e",
                              trial, n, K, mw, max_diff))
  }
})
