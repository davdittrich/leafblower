test_that("A7: all 5 quality metrics present in calib_result for oris", {
  set.seed(42)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "oris",
                           max_iterations = 500,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error", "mean_error", "kl", "chi2", "l1_weight_change"))
    expect_true(nm %in% names(result),
                info = sprintf("metric '%s' missing from calib_result", nm))
  expect_true(is.finite(result$max_error))
  expect_true(is.finite(result$mean_error))
  expect_true(is.finite(result$kl))
  expect_true(is.finite(result$chi2) || is.infinite(result$chi2))  # chi2 can be Inf on degenerate
  expect_true(is.finite(result$l1_weight_change))
})

test_that("A7: all 5 quality metrics present in calib_result for raking", {
  set.seed(43)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3"), n, replace = TRUE))
  )
  target <- list(a = c("1"=0.5,"2"=0.5), b = c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "raking",
                           max_iterations = 500,
                           convergence = list(absolute = 1e-6),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  for (nm in c("max_error", "mean_error", "kl", "chi2", "l1_weight_change"))
    expect_true(nm %in% names(result))
  expect_true(is.finite(result$mean_error))
  expect_true(is.finite(result$l1_weight_change))
})


test_that("A7: metrics non-zero after max_iter exit (solver exits before kErrCheckInterval)", {
  # Use max_iterations=1 to force exit after 1 iteration, likely before kErrCheckInterval.
  # All three solvers check at iter==1, so metrics must be populated.
  set.seed(55)
  n <- 500
  data <- data.frame(a = factor(sample(c("1","2","3"), n, replace = TRUE)))
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- leafblower::harvest(data, target, max_weight = 5, method = "oris",
                           max_iterations = 1,
                           convergence = list(absolute = 1e-20),  # impossible threshold
                           attach_weights = FALSE)
  result <- attr(w, "result")
  # Metrics must be populated (finite) even after 1 iteration:
  expect_true(is.finite(result$mean_error))
  expect_true(is.finite(result$l1_weight_change))
  expect_true(result$l1_weight_change >= 0)
})

test_that("a0gk: metrics finite at exit with MAX_ERR criterion (gated path)", {
  # Gate: mean_err/kl/chi2 skipped at intermediate checks when criterion=MAX_ERR/PCT.
  # They are computed on the check where convergence fires (about_to_converge gate)
  # and on the final budget iteration. Verify all three metrics are finite at exit.
  # Use 2 imbalanced margins so calibration takes multiple iterations.
  set.seed(7)
  n <- 600
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE, prob = c(0.6, 0.3, 0.1))),
    b = factor(sample(c("X","Y"),     n, replace = TRUE, prob = c(0.4, 0.6)))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("X"=0.5,"Y"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "oris",
                           max_iterations = 100,
                           convergence = list(absolute = 1e-3),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_true(is.finite(result$mean_error))
  expect_true(is.finite(result$kl) && result$kl >= 0)
  expect_true(is.finite(result$chi2) && result$chi2 >= 0)
  # Convergence must have fired (not budget exhaustion) to exercise the gate path
  expect_true(result$iterations < 100)
})

test_that("B2 dispatch: NA data uses K-pass (not single-pass), no crash", {
  # K=3 with NA in margin col → !anyNA fails → fallback to K-pass.
  set.seed(42L); n <- 500L
  df_na <- data.frame(
    a = factor(sample(c("x","y","z"), n, TRUE)),
    b = factor(sample(c("M","F"), n, TRUE)),
    c = factor(sample(c("1","2","3"), n, TRUE))
  )
  df_na$a[sample(n, 25L)] <- NA
  tgt_na <- list(
    a = c(x=0.5, y=0.3, z=0.2),
    b = c(M=0.6, F=0.4),
    c = c("1"=0.4, "2"=0.4, "3"=0.2)
  )
  w_na <- rep(1, n)
  qm <- leafblower:::compute_quality_metrics(w_na, tgt_na, df_na)
  # Dispatch to K-pass must succeed without error and return finite/Inf
  expect_true(is.finite(qm$margin_kl) || is.infinite(qm$margin_kl))
  expect_true(is.finite(qm$design_effect))
  expect_true(is.finite(qm$weight_kl))
})

test_that("B2 values: single-pass equals K-pass on no-NA K>=3 data", {
  # K=3 no-NA → single-pass path. Verify numerical equivalence vs K-pass.
  set.seed(42L); n <- 500L
  df_clean <- data.frame(
    a = factor(sample(c("x","y","z"), n, TRUE)),
    b = factor(sample(c("M","F"), n, TRUE)),
    c = factor(sample(c("1","2","3"), n, TRUE))
  )
  tgt_c <- list(
    a = c(x=0.5, y=0.3, z=0.2),
    b = c(M=0.6, F=0.4),
    c = c("1"=0.4, "2"=0.4, "3"=0.2)
  )
  w_c <- runif(n, 0.5, 2.0)
  qm <- leafblower:::compute_quality_metrics(w_c, tgt_c, df_clean)

  # Reference: K-pass computation inline (mirrors fallback branch)
  Z <- sum(w_c)
  ref_margin_kl <- sum(sapply(names(tgt_c), function(k) {
    T_k <- tgt_c[[k]]
    obs_k <- df_clean[[k]]
    valid <- !is.na(obs_k)
    Z_k <- sum(w_c[valid])
    W_k <- tapply(w_c[valid], droplevels(obs_k[valid]), sum) / Z_k
    if (any(T_k[setdiff(names(T_k), names(W_k))] > 0)) return(Inf)
    common <- intersect(names(T_k), names(W_k))
    T_sub <- T_k[common]; W_sub <- W_k[common]
    pos <- T_sub > 0
    if (!any(pos)) return(0)
    sum(T_sub[pos] * log(T_sub[pos] / pmax(W_sub[pos], 1e-15)))
  }))
  expect_equal(qm$margin_kl, ref_margin_kl, tolerance = 1e-12)
  expect_true(is.finite(qm$margin_kl))
})

test_that("B2 dispatch: K<3 uses K-pass (single-pass overhead unprofitable)", {
  # K=2 → length(target_list) >= 3 fails → fallback to K-pass.
  set.seed(43L); n <- 300L
  df_k2 <- data.frame(
    a = factor(sample(c("x","y","z"), n, TRUE)),
    b = factor(sample(c("M","F"), n, TRUE))
  )
  tgt_k2 <- list(a = c(x=0.4, y=0.35, z=0.25), b = c(M=0.5, F=0.5))
  w_k2 <- rep(1, n)
  qm <- leafblower:::compute_quality_metrics(w_k2, tgt_k2, df_k2)
  expect_true(is.finite(qm$margin_kl) && qm$margin_kl >= 0)
})

test_that("B1: compute_quality_metrics extraction: values identical to inline", {
  # Snapshot test capturing expected values from inline block (lines 536-575)
  # before extraction to helper. After extraction, verify helper produces identical results.
  set.seed(11L)
  n <- 1000L
  df_t <- data.frame(
    a = factor(sample(c("x","y","z"), n, TRUE, prob = c(0.5, 0.3, 0.2))),
    b = factor(sample(c("M","F"), n, TRUE, prob = c(0.55, 0.45)))
  )
  tgt_t <- list(a = c(x=0.4, y=0.35, z=0.25), b = c(M=0.5, F=0.5))
  r_base <- leafblower::harvest(df_t, tgt_t, method = "oris", max_weight = 5,
                                max_iterations = 50L, attach_weights = FALSE)
  w <- r_base
  res <- attr(r_base, "result")
  expected_margin_kl <- res$margin_kl
  expected_deff <- res$design_effect
  expected_weight_kl <- res$weight_kl
  expected_eff_obs <- res$effective_observations

  # leafblower-on7a v3 boundary: asserts qm$design_effect (Kish) == res$design_effect (Kish).
  # Neither invokes the H&V 4-arg design_effect(); that is intentionally not plumbed here.
  # Call extracted helper and verify identical values
  qm <- leafblower:::compute_quality_metrics(w, tgt_t, df_t)
  expect_equal(qm$margin_kl, expected_margin_kl, tolerance = 1e-10)
  expect_equal(qm$design_effect, expected_deff, tolerance = 1e-10)
  expect_equal(qm$weight_kl, expected_weight_kl, tolerance = 1e-10)
  expect_equal(qm$effective_observations, expected_eff_obs, tolerance = 1e-10)
})

test_that("mb06: compute_quality_metrics margin_kl finite for CHARACTER margin; bit-equals factor path", {
  # Regression (function under test = compute_quality_metrics, NOT the solver):
  # droplevels() has no method for class character -> the K-pass fallback errored
  # ("no applicable method for 'droplevels' ...") -> tryCatch swallowed it ->
  # margin_kl = NA + warning. A CHARACTER margin forces the K-pass path
  # (single-pass requires ALL margin cols be factors). The fix guards droplevels()
  # to factors; tapply coerces the character index to a factor over present values
  # only -> bit-identical W_k -> bit-identical margin_kl. Call the helper directly
  # with a fixed weights vector so NO solver warning can be misattributed here.
  set.seed(2026L)
  n     <- 800L
  cats  <- sample(c("x", "y", "z"), n, TRUE, prob = c(0.5, 0.3, 0.2))
  w     <- runif(n, 0.5, 2.0)                                # fixed positive weights
  tgt   <- list(a = c(x = 0.4, y = 0.35, z = 0.25))
  df_chr <- data.frame(a = cats, stringsAsFactors = FALSE)   # CHARACTER margin -> K-pass
  df_fac <- data.frame(a = factor(cats))                     # same data, FACTOR margin

  # Only the helper call is wrapped: no solver runs, so no solver warning is possible.
  qm_chr <- expect_no_warning(leafblower:::compute_quality_metrics(w, tgt, df_chr))
  mk_chr <- qm_chr$margin_kl
  expect_true(is.finite(mk_chr))                             # (a) finite / NOT NA

  qm_fac <- leafblower:::compute_quality_metrics(w, tgt, df_fac)
  mk_fac <- qm_fac$margin_kl
  # (c) char path feeds bit-identical inputs to the KL helper as the factor path
  expect_identical(mk_chr, mk_fac)
})

test_that("dtkn.12: add_na_proportion 'NA' bin -> FINITE margin_kl (NA obs counted)", {
  # Regression (function under test = compute_quality_metrics, NOT the solver).
  # When the target carries an injected "NA" level (add_na_proportion), the
  # K-pass fallback used to exclude is.na(obs) -> W_k had no "NA" level while
  # T_k did -> margin_kl_one() returned +Inf. Fix: recode NA -> "NA" and
  # normalize over ALL obs so W_k gains an "NA" mass ~ target's NA fraction.
  set.seed(7L)
  n   <- 400L
  a   <- sample(c("x", "y", "z"), n, TRUE, prob = c(0.5, 0.3, 0.2))
  a[sample(n, 40L)] <- NA                              # 10% NA
  df  <- data.frame(a = factor(a))                     # NA kept as NA
  w   <- runif(n, 0.5, 2.0)                            # fixed positive weights
  tgt <- list(a = c(x = 0.45, y = 0.27, z = 0.18, "NA" = 0.10))  # sums to 1

  qm <- leafblower:::compute_quality_metrics(w, tgt, df)
  expect_true(is.finite(qm$margin_kl))                 # (a) finite, was +Inf

  # (b) NA-bin contributes ~0 to KL: weighted NA mass ≈ target 0.10, so the
  # NA term T*log(T/W) is small. Full margin_kl stays modest, not blown up.
  expect_lt(qm$margin_kl, 1.0)
})

test_that("dtkn.12: no 'NA' target level -> NA obs still EXCLUDED (unchanged path)", {
  # Guard: the fix keys ONLY on a target level literally named "NA". A margin
  # with NA data but NO "NA" target level must keep the old exclude-NA behavior
  # (normalize over non-NA obs), bit-identical to pre-fix.
  set.seed(11L)
  n   <- 300L
  a   <- sample(c("x", "y", "z"), n, TRUE, prob = c(0.5, 0.3, 0.2))
  a[sample(n, 30L)] <- NA
  df  <- data.frame(a = factor(a))
  w   <- runif(n, 0.5, 2.0)
  tgt <- list(a = c(x = 0.5, y = 0.3, z = 0.2))        # NO "NA" level

  qm <- leafblower:::compute_quality_metrics(w, tgt, df)
  # Reconstruct the exclude-NA expectation independently:
  valid <- !is.na(a); Z <- sum(w[valid])
  Wk <- tapply(w[valid], droplevels(factor(a[valid])), sum) / Z
  T_k <- tgt$a
  exp_kl <- sum(T_k * log(T_k / Wk[names(T_k)]))
  expect_equal(qm$margin_kl, exp_kl, tolerance = 1e-12)
})
