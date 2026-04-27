test_that("T1a: chebyshev and grake are now implemented (no stub error)", {
  data <- data.frame(a = factor(c("1","2","1","2","1")))
  target <- list(a = c("1"=0.5, "2"=0.5))
  expect_no_error(
    leafblower::harvest(data, target, max_weight=3, method="chebyshev", attach_weights=FALSE)
  )
  expect_no_error(
    leafblower::harvest(data, target, max_weight=3, method="grake", attach_weights=FALSE)
  )
})

test_that("T1b: convergence_used$objective and $minimized_metric present", {
  set.seed(1)
  data <- data.frame(a=factor(sample(c("1","2"),200,TRUE)))
  target <- list(a=c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=3, method="ieppa",
                           attach_weights=FALSE)
  r <- attr(w, "result")
  expect_true("objective" %in% names(r$convergence_used))
  expect_true("minimized_metric" %in% names(r$convergence_used))
  expect_true(is.finite(r$convergence_used$objective))
})

test_that("T3: ieppa default convergence is kl+improvement", {
  set.seed(42)
  data <- data.frame(
    a = factor(sample(c("1","2","3"), 500, replace=TRUE)),
    b = factor(sample(c("1","2"), 500, replace=TRUE))
  )
  target <- list(a=c("1"=1/3,"2"=1/3,"3"=1/3), b=c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=10, method="ieppa",
                           max_iterations=500, attach_weights=FALSE)
  r <- attr(w, "result")
  expect_equal(r$convergence_used$metric, "kl",
               info="ieppa default metric must be 'kl' after spec change")
  expect_equal(r$convergence_used$rule, "improvement")
})

test_that("A1: sinkhorn KL <= ieppa KL at best_iter", {
  skip_on_cran()
  skip_if_not_installed("arrow")
  skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  ref_path <- test_path("fixtures/ieppa_kl_reference_stepstone.rds")
  skip_if(!file.exists(ref_path), "ieppa KL reference fixture not generated yet")
  ref <- readRDS(ref_path)

  # No unconditional skip — test FAILS with "sinkhorn not yet implemented"
  # until Plan C implements method="sinkhorn". That failure IS the RED state.
  data <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  tgt_raw <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
  target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })

  # Plan A: this expect_no_error call FAILS (sinkhorn is a stub returning RK_ERR_BADARG)
  # Plan C GREEN: sinkhorn is implemented so expect_no_error passes
  expect_no_error(
    w_s <- leafblower::harvest(data, target, method="sinkhorn",
                                max_weight=5, max_iterations=3000,
                                attach_weights=FALSE),
    message="Plan C must implement method='sinkhorn' to reach this point"
  )
  r_s <- attr(w_s, "result")
  expect_lte(r_s$convergence_used$objective, ref$kl_at_best_iter,
             label="sinkhorn KL <= ieppa best_iter KL")
  expect_equal(r_s$status, 0L, label="sinkhorn must converge")
  expect_equal(r_s$convergence_used$fired_at_iter, r_s$iterations,
               label="A7: sinkhorn best_iter == last_iter (monotone)")
})

test_that("T-routing: AUTO uses raking when incompressible (M_cell/n > 0.9)", {
  # 200 obs each in its own group → M_cell = 200, ratio = 1.0 (fully incompressible).
  # AUTO must select raking (not iEPPA) because there is no cell compression benefit.
  set.seed(99)
  n2 <- 200
  # Each obs is its own group for margin 'a' → M_cell = n2
  data <- data.frame(
    a = factor(seq_len(n2)),
    b = factor(sample(c("1", "2"), n2, replace = TRUE))
  )
  target <- list(
    a = setNames(rep(1 / n2, n2), as.character(seq_len(n2))),
    b = c("1" = 0.5, "2" = 0.5)
  )
  w <- leafblower::harvest(data, target, max_weight = 5, method = "auto",
                           max_iterations = 200, attach_weights = FALSE)
  r <- attr(w, "result")
  expect_equal(r$status, 0L,
               label = "AUTO must converge even with incompressible data")
  expect_lt(r$max_error, 0.1)
  # Routing check: RK_ALG_RAKING == 3
  expect_equal(r$algorithm_used, 3L,
               label = "AUTO must select raking when M_cell/n > 0.9")
})

test_that("A5: raking cell-table: correct + speedup vs obs-level reference", {
  skip_on_cran()
  skip_if_not_installed("arrow")
  skip_if_not_installed("jsonlite")
  skip_if(!file.exists("benchmarks/stepstone_fulldata_bench_data.parquet"),
          "stepstone benchmark data not available")
  ref_path <- test_path("fixtures/raking_obs_reference_stepstone.rds")
  skip_if(!file.exists(ref_path), "obs-level raking reference fixture not generated")
  ref <- readRDS(ref_path)

  library(arrow); library(jsonlite)
  data <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
  tgt_raw <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
  target <- lapply(tgt_raw, function(x) { v <- unlist(x); v / sum(v) })

  t0 <- proc.time()["elapsed"]
  suppressWarnings(
    w_cell <- leafblower::harvest(data, target, method="raking",
                                  max_weight=5, max_iterations=500, attach_weights=FALSE)
  )
  elapsed_cell <- proc.time()["elapsed"] - t0

  r_cell <- attr(w_cell, "result")

  # Convergence
  expect_equal(r_cell$status, 0L, label="raking cell-table must converge")

  # Correctness: max_err within 1e-8 of obs-level reference (homogeneous d_i)
  expect_lt(abs(r_cell$max_error - ref$max_error), 1e-8,
            label="cell-table max_err must match obs-level to 1e-8")

  # Speedup: < 5% of obs-level elapsed (spec A5: elapsed_cell < elapsed_obs * 0.05)
  expect_lt(elapsed_cell, ref$elapsed_obs * 0.05,
            label=sprintf("cell-table must be < 5%% of obs-level elapsed (obs=%.1fs)", ref$elapsed_obs))

  # Hard bounds: all weights in [min_weight, max_weight]
  expect_true(all(w_cell >= 1/5 - 1e-10 & w_cell <= 5 + 1e-10),
              label="all weights within [1/5, 5] bounds")
})

test_that("R-bounds: raking respects min_weight/max_weight exactly after fix", {
  # Tight bounds + skewed targets force clamping; normalization after clamp would violate bounds.
  set.seed(7)
  n <- 1000
  data <- data.frame(
    a = factor(sample(c("1","2","3","4","5"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
  target <- list(
    a = c("1"=0.5, "2"=0.2, "3"=0.15, "4"=0.1, "5"=0.05),
    b = c("1"=0.8, "2"=0.2)
  )
  out <- leafblower::harvest(data, target, min_weight=0.5, max_weight=1.5,
                             method="raking", max_iterations=500, attach_weights=TRUE)
  w <- out$weights
  expect_true(all(w >= 0.5 - 1e-10),
              info=sprintf("min weight %.6f < 0.5", min(w)))
  expect_true(all(w <= 1.5 + 1e-10),
              info=sprintf("max weight %.6f > 1.5", max(w)))
})

test_that("T-overflow: AUTO routing + algorithm_used populated", {
  set.seed(1)
  n <- 1e5L
  data <- data.frame(
    a = factor(sample(1:100, n, replace=TRUE)),
    b = factor(sample(1:2, n, replace=TRUE))
  )
  target <- list(
    a = setNames(rep(0.01, 100), as.character(1:100)),
    b = c("1"=0.5, "2"=0.5)
  )
  w <- leafblower::harvest(data, target, max_weight=5, method="auto", attach_weights=FALSE)
  r <- attr(w, "result")
  expect_equal(r$status, 0L)
  expect_lt(r$max_error, 0.01)
  # algorithm_used: AUTO with M_cell << n should route to ieppa (algorithm=1)
  expect_equal(r$algorithm_used, 1L,
               info="AUTO must select ieppa when M_cell/n << 0.9")
})

test_that("T-auto-kl: method='auto' defaults to kl convergence metric", {
  set.seed(3)
  data <- data.frame(a=factor(sample(c("1","2","3"),300,TRUE)))
  target <- list(a=c("1"=0.4,"2"=0.4,"3"=0.2))
  w <- leafblower::harvest(data, target, max_weight=5, method="auto", attach_weights=FALSE)
  r <- attr(w, "result")
  expect_equal(r$convergence_used$metric, "kl",
               info="AUTO must use kl default — both ieppa and raking default to kl")
})

test_that("S1: sinkhorn handles tight bounds without overflow (a[c] clamp guard)", {
  # Regression guard: without a[c] clamping, exp overflows after ~700 iters of persistent clamping.
  set.seed(99)
  n <- 500
  data <- data.frame(
    a = factor(sample(c("1","2","3","4","5"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
  target <- list(
    a = c("1"=0.6, "2"=0.2, "3"=0.1, "4"=0.05, "5"=0.05),
    b = c("1"=0.9, "2"=0.1)
  )
  w <- leafblower::harvest(data, target, min_weight=0.1, max_weight=2.0,
                           method="sinkhorn", max_iterations=2000, attach_weights=FALSE)
  r <- attr(w, "result")
  expect_false(r$status == 3L,
               info=sprintf("sinkhorn returned INFEAS (status=%d) — likely a[c] overflow", r$status))
  expect_true(r$status %in% c(0L, 1L),
              info=sprintf("expected 0 (OK) or 1 (NOCONV), got %d", r$status))
  expect_true(all(w >= 0.1 - 1e-10 & w <= 2.0 + 1e-10),
              info="bounds violated")
})

test_that("S2: sinkhorn achieves KL <= raking on synthetic (no external data)", {
  set.seed(7)
  n <- 300
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
  target <- list(
    a = c("1"=0.5, "2"=0.3, "3"=0.2),
    b = c("1"=0.6, "2"=0.4)
  )
  w_sink <- leafblower::harvest(data, target, max_weight=5, method="sinkhorn",
                                max_iterations=500, attach_weights=FALSE)
  w_rake <- leafblower::harvest(data, target, max_weight=5, method="raking",
                                max_iterations=500, attach_weights=FALSE)
  r_sink <- attr(w_sink, "result")
  r_rake <- attr(w_rake, "result")
  expect_equal(r_sink$status, 0L, info="sinkhorn must converge")
  expect_equal(r_rake$status, 0L, info="raking must converge")
  expect_lte(r_sink$kl, r_rake$kl + 1e-6,
             label="sinkhorn KL <= raking KL (competitive KL minimizer)")
  expect_true(all(w_sink >= 1/5 - 1e-10 & w_sink <= 5 + 1e-10),
              info="sinkhorn weights must respect bounds")
})

test_that("D1: greg achieves chi2 <= raking on synthetic (no external data)", {
  set.seed(5)
  n <- 400
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
  target <- list(a=c("1"=0.4,"2"=0.4,"3"=0.2), b=c("1"=0.6,"2"=0.4))
  w_greg <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                                method="greg", attach_weights=FALSE)
  w_rake <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                                method="raking", max_iterations=500, attach_weights=FALSE)
  r_greg <- attr(w_greg, "result")
  r_rake <- attr(w_rake, "result")
  expect_equal(r_greg$status, 0L, info="greg must converge")
  expect_equal(r_rake$status, 0L, info="raking must converge")
  expect_lte(r_greg$chi2, r_rake$chi2 + 1e-6,
             label="greg chi2 <= raking chi2")
  expect_true(all(w_greg >= 0.2 - 1e-10 & w_greg <= 5 + 1e-10),
              info="greg bounds must hold")
})

test_that("E1: chebyshev max_err <= raking max_err (correctness)", {
  set.seed(11)
  n <- 400
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
  target <- list(a=c("1"=0.4,"2"=0.4,"3"=0.2), b=c("1"=0.6,"2"=0.4))
  w_cheb <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                                method="chebyshev", attach_weights=FALSE)
  r_cheb <- attr(w_cheb, "result")
  w_rake <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                                method="raking", max_iterations=500, attach_weights=FALSE)
  r_rake <- attr(w_rake, "result")
  expect_equal(r_cheb$status, 0L)
  expect_lte(r_cheb$max_error, r_rake$max_error + 1e-6,
             label="chebyshev max_err <= raking max_err")
  expect_true(all(w_cheb >= 0.2 - 1e-10 & w_cheb <= 5 + 1e-10))
})

test_that("E2: grake grake_norm <= raking grake_norm (correctness)", {
  set.seed(13)
  n <- 400
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace=TRUE)),
    b = factor(sample(c("1","2"), n, replace=TRUE))
  )
  target <- list(a=c("1"=0.4,"2"=0.4,"3"=0.2), b=c("1"=0.6,"2"=0.4))
  w_grake <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                                 method="grake", attach_weights=FALSE)
  r_grake <- attr(w_grake, "result")
  w_rake <- leafblower::harvest(data, target, min_weight=0.2, max_weight=5,
                                method="raking", max_iterations=500, attach_weights=FALSE)
  r_rake <- attr(w_rake, "result")
  expect_equal(r_grake$status, 0L)
  expect_lte(r_grake$grake_norm, r_rake$grake_norm + 1e-6,
             label="grake grake_norm <= raking grake_norm")
  expect_true(all(w_grake >= 0.2 - 1e-10 & w_grake <= 5 + 1e-10))
})

# ──────────────────────────────────────────────────────────────────────────────
# T1.B regression: no linear overflow on skewed multi-margin problems.
# Uses K=20, n=100000 (large enough to avoid S_lin collapse), max_weight=3.
# Math: kLinearOverflowTrip = (DBL_MAX/2)^(1/20) ≈ 2.1e15.
# Product across K=20 margins grows as (1.5)^(20*N) per sweep.
# Before T1.B: "overflow trip" message fires at ~iter 4-5.
# After  T1.B: renorm fires, no overflow message.
# Note: max_weight=3 prevents convergence in 50 iters (tight bounds + oscillation)
# so we test absence of overflow message, NOT status == 0.
# ──────────────────────────────────────────────────────────────────────────────
test_that("ieppa: no linear overflow trip on K=20 skewed targets (T1.B)", {
  set.seed(42); K <- 20L; n <- 100000L
  cols <- paste0("v", seq_len(K))
  data <- as.data.frame(lapply(seq_len(K),
    function(k) factor(sample(5L, n, TRUE), levels = 1:5)))
  names(data) <- cols
  skewed <- c("1" = 0.3, "2" = 0.175, "3" = 0.175, "4" = 0.175, "5" = 0.175)
  target <- setNames(lapply(seq_len(K), function(.) skewed), cols)

  # Capture C++ stdout where st.log() writes (overflow message goes here)
  out <- capture.output(
    r <- leafblower::harvest(data, target, method = "ieppa",
                             min_weight = 0.2, max_weight = 3,
                             max_iterations = 50,
                             attach_weights = FALSE, verbose = 1),
    type = "output"
  )
  overflow_msgs <- grep("overflow trip", out, value = TRUE)

  # Before T1.B: "overflow trip" message fires → overflow_msgs length > 0 → FAIL
  # After  T1.B: renorm prevents overflow → no message → PASS
  expect_length(overflow_msgs, 0L)
})
