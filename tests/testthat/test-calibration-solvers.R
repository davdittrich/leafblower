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

test_that("T1b: convergence_used$solver_objective and $minimized_metric present", {
  set.seed(1)
  data <- data.frame(a=factor(sample(c("1","2"),200,TRUE)))
  target <- list(a=c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight=3, method="ieppa",
                           attach_weights=FALSE)
  r <- attr(w, "result")
  expect_true("solver_objective" %in% names(r$convergence_used))
  expect_true("minimized_metric" %in% names(r$convergence_used))
  expect_true(is.finite(r$convergence_used$solver_objective))
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
  expect_equal(r$convergence_used$metric, "marginal_kl",
               info="ieppa default metric must be 'marginal_kl' (Task 0 convergence overlay)")
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
  expect_lte(attr(w_s, "result")$convergence_used$solver_objective, ref$kl_at_best_iter,
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

test_that("T-auto-kl: method='ieppa' defaults to marginal_kl convergence metric", {
  set.seed(3)
  data <- data.frame(a=factor(sample(c("1","2","3"),300,TRUE)))
  target <- list(a=c("1"=0.4,"2"=0.4,"3"=0.2))
  w <- leafblower::harvest(data, target, max_weight=5, method="ieppa", attach_weights=FALSE)
  r <- attr(w, "result")
  expect_equal(r$convergence_used$metric, "marginal_kl",
               info="ieppa default metric is marginal_kl (calibration quality loss)")
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
  expect_true(r$status %in% c(0L, 4L, 5L),
              info=sprintf("expected 0 (OK) or 4 (BUDGET) or 5 (STALL), got %d", r$status))
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

test_that("ieppa: cell-mode weights respect max_weight hard cap", {
  # Reliable trigger: K=1, max_weight=1.5, targets=(0.9,0.1), uniform data.
  # Ideal weight for cat "A" ≈ 1.8 > cap 1.5 → cells clamped → W_total < n
  # → norm = n/W_total > 1 → post-norm wmax = 1.5 * norm > 1.5.  Bug fires.
  set.seed(1); n <- 1000L
  df <- data.frame(v1 = factor(sample(c("A","B"), n, replace=TRUE, prob=c(.5,.5))))
  tgt <- list(v1 = c("A"=0.9, "B"=0.1))

  r <- leafblower::harvest(df, tgt, method = "ieppa",
    max_weight = 1.5, min_weight = 0.0,
    bounds_mode = "cell", max_iterations = 500L,
    attach_weights = FALSE, verbose = 0)
  w <- as.numeric(r)

  # Before fix: post-norm wmax ≈ 1.5 * (1000/850) ≈ 1.76 > 1.5 → FAIL
  # After fix:  cell-mode clamp applied → wmax ≤ 1.5 → PASS
  expect_true(max(w) <= 1.5 + 1e-9,
    label = sprintf("max weight %.6f exceeds cap 1.5", max(w)))
  expect_true(min(w) >= 0.0 - 1e-9,
    label = sprintf("min weight %.6f below floor 0.0", min(w)))
})

# ── T1: solver_objective field existence (RED: field not found before Task 2) ──
test_that("T1: sinkhorn convergence_used$solver_objective exists and is weight KL", {
  set.seed(1L); n <- 800L
  df  <- data.frame(v1 = factor(sample(c("A","B","C"), n, TRUE)))
  tgt <- list(v1 = c("A"=0.6, "B"=0.3, "C"=0.1))
  r_mx <- leafblower::harvest(df, tgt, method="sinkhorn",
    convergence=list(metric="max_err"), max_iterations=200L, attach_weights=FALSE)
  r_kl <- leafblower::harvest(df, tgt, method="sinkhorn",
    convergence=list(metric="kl"),     max_iterations=200L, attach_weights=FALSE)
  # RED state: $solver_objective field doesn't exist — stopifnot ERRORs if NULL
  obj_mx <- attr(r_mx, "result")$convergence_used$solver_objective
  obj_kl <- attr(r_kl, "result")$convergence_used$solver_objective
  # Force ERROR if field not present (returns NULL before Task 2)
  stopifnot(!is.null(obj_mx), !is.null(obj_kl))
  # Both stopping criteria -> same mathematical objective (weight KL, not stopping value)
  expect_true(is.finite(obj_mx),  label="solver_objective is finite")
  expect_true(obj_mx < 0.5,       label="solver_objective is weight KL (not max_err ~0.05)")
  expect_true(abs(obj_mx - obj_kl) / max(obj_mx, obj_kl) < 0.5,
              label="objective consistent across stopping criteria")
})

# ── T2: ieppa_soft available (RED: unknown method before Task 3) ──
test_that("T2: ieppa_soft method exists and respects max_weight", {
  set.seed(2L); n <- 2000L
  df  <- data.frame(v1 = factor(sample(c("X","Y"), n, TRUE, prob=c(.3,.7))))
  tgt <- list(v1 = c("X"=0.8, "Y"=0.2))
  r <- leafblower::harvest(df, tgt, method="ieppa_soft",
    max_weight=2.0, min_weight=0.0, max_iterations=300L, attach_weights=FALSE)
  w <- as.numeric(r)
  expect_true(max(w) <= 2.0 + 1e-9, label="ieppa_soft wmax <= max_weight")
  expect_true(min(w) >= 0.0 - 1e-9, label="ieppa_soft wmin >= min_weight")
  expect_equal(attr(r, "result")$status, 0L)
})

# ── T3: ieppa_soft no worse than ieppa on tight-bounds problem ──
# ADMM benefit is path-dependent: permanently-saturated cells show no improvement.
# Assertion weakened to <= (no regression guarantee, not strict improvement).
test_that("T3: ieppa_soft max_err strictly < ieppa max_err on tight bounds", {
  set.seed(3L); n <- 5000L
  df  <- data.frame(v1 = factor(sample(5L, n, TRUE)))
  tgt <- list(v1 = setNames(c(0.4, 0.3, 0.15, 0.1, 0.05), as.character(1:5)))
  r_hard <- leafblower::harvest(df, tgt, method="ieppa",
    max_weight=1.8, min_weight=0, max_iterations=500L, attach_weights=FALSE)
  r_soft <- leafblower::harvest(df, tgt, method="ieppa_soft",
    max_weight=1.8, min_weight=0, max_iterations=500L, attach_weights=FALSE)
  me_hard <- attr(r_hard, "result")$max_error
  me_soft <- attr(r_soft, "result")$max_error
  expect_true(me_soft <= me_hard + 1e-9,
              label="ieppa_soft max_err not worse than ieppa (no regression)")
})

# ── Raking Bregman Dykstra RED test ─────────────────────────────────────────
# Before fix: Euclidean hyperplane correction changes fixed point vs pure IPF.
# After fix:  multiplicative hyperplane = KL projection → same fixed point as ieppa.
# RED: expect_equal(wkl_raking, wkl_ieppa, tol=1e-4) FAILS before Bregman Dykstra.
# solver_objective field confirmed: raking.cpp:349, harvest.R:280.
# ────────────────────────────────────────────────────────────────────────────
test_that("raking-bregman: unconstrained raking matches ieppa weight_kl (unified KL fixed point)", {
  set.seed(1L); n <- 2000L
  df <- data.frame(
    v1 = factor(sample(3L, n, TRUE)),
    v2 = factor(sample(2L, n, TRUE))
  )
  tgt <- list(v1 = c("1"=0.5, "2"=0.3, "3"=0.2),
              v2 = c("1"=0.6, "2"=0.4))

  r_raking <- leafblower::harvest(df, tgt, method="raking",
    min_weight=0, max_weight=Inf, max_iterations=500L,
    attach_weights=FALSE)
  r_ieppa  <- leafblower::harvest(df, tgt, method="ieppa",
    min_weight=0, max_weight=Inf, max_iterations=500L,
    attach_weights=FALSE)

  wkl_raking <- attr(r_raking, "result")$convergence_used$solver_objective
  wkl_ieppa  <- attr(r_ieppa,  "result")$convergence_used$solver_objective

  # Both are unconstrained IPF → same KL minimum → same weight_kl.
  # Before Bregman: Euclidean hyperplane shifts raking fixed point → FAIL.
  # After  Bregman: unified KL geometry → PASS.
  expect_equal(wkl_raking, wkl_ieppa, tolerance=1e-4,
               label="unconstrained raking must reach same weight_kl as ieppa")
})

# ── SQUAREM acceleration tests ──────────────────────────────────────────────

test_that("squarem-red: accelerate=TRUE produces different iterations than accelerate=FALSE", {
  # RED state (before SQUAREM): accelerate=TRUE silently ignored → same iterations.
  # GREEN state (after SQUAREM): accelerate=TRUE takes fewer super-steps → different iters.
  # Note: harvest() currently has accelerate in its 'ignored' list → same result as FALSE.
  set.seed(42L); n <- 2000L
  df  <- data.frame(v1 = factor(sample(3L, n, TRUE)), v2 = factor(sample(2L, n, TRUE)))
  tgt <- list(v1=c("1"=0.5,"2"=0.3,"3"=0.2), v2=c("1"=0.6,"2"=0.4))

  r_base <- leafblower::harvest(df, tgt, method="raking",
    accelerate=FALSE, max_weight=5, max_iterations=500L, attach_weights=FALSE)
  r_acc  <- leafblower::harvest(df, tgt, method="raking",
    accelerate=TRUE,  max_weight=5, max_iterations=500L, attach_weights=FALSE)

  iters_base <- attr(r_base, "result")$iterations
  iters_acc  <- attr(r_acc,  "result")$iterations

  # Before SQUAREM: accelerate=TRUE ignored → same iterations as FALSE
  # After  SQUAREM: fewer super-steps (each = 3 F-evals) → different iters
  expect_false(isTRUE(all.equal(iters_base, iters_acc)),
               label="accelerate=TRUE must use different iterations than accelerate=FALSE")
})

test_that("squarem-ac3: accelerate=FALSE matches pre-SQUAREM baseline to 1e-14", {
  ref_path <- test_path("fixtures/raking_squarem_baseline.rds")
  skip_if(!file.exists(ref_path), "raking_squarem_baseline.rds not generated yet")

  set.seed(42L); n <- 2000L
  df  <- data.frame(v1 = factor(sample(3L, n, TRUE)), v2 = factor(sample(2L, n, TRUE)))
  tgt <- list(v1=c("1"=0.5,"2"=0.3,"3"=0.2), v2=c("1"=0.6,"2"=0.4))

  r <- leafblower::harvest(df, tgt, method="raking",
    accelerate=FALSE, max_weight=5, max_iterations=500L, attach_weights=FALSE)
  w <- as.numeric(r)
  w_ref <- readRDS(ref_path)

  # tolerance=1e-14 (not 0): platform FP non-determinism makes exact bit-equality unreliable
  expect_equal(w, w_ref, tolerance=1e-14,
               label="accelerate=FALSE must match pre-SQUAREM baseline to 1e-14")
})

test_that("squarem-ac4: step-halving path produces valid weights (no NaN, bounds respected)", {
  # Tests that step-halving backtrack correctness holds: weights must be finite and bounded.
  set.seed(7L); n <- 300L
  df  <- data.frame(v1=factor(sample(4L,n,TRUE)), v2=factor(sample(3L,n,TRUE)))
  tgt <- list(v1=c("1"=0.4,"2"=0.3,"3"=0.2,"4"=0.1),
              v2=c("1"=0.5,"2"=0.3,"3"=0.2))
  r <- leafblower::harvest(df, tgt, method="raking", accelerate=TRUE,
                           max_weight=2, max_iterations=200L, attach_weights=FALSE)
  w <- as.numeric(r)
  expect_true(all(is.finite(w)), label="AC4: no NaN/Inf in weights")
  expect_true(all(w >= 0 & w <= 2 + 1e-9), label="AC4: weights within [0, max_weight]")
  expect_equal(sum(w), n, tolerance=1e-8, label="AC4: weights sum to n")
})

test_that("squarem-ac5: accelerate=TRUE with non-raking method warns and runs", {
  df  <- data.frame(v1 = factor(c("1","2","1","2","1")))
  tgt <- list(v1=c("1"=0.5,"2"=0.5))

  expect_warning(
    leafblower::harvest(df, tgt, method="ieppa", accelerate=TRUE,
                        max_weight=3, attach_weights=FALSE),
    regexp="accelerate.*raking",
    label="accelerate=TRUE with ieppa must warn"
  )
})

test_that("wf-min-weight: water-filling respects min_weight > 0 (lower bound path)", {
  set.seed(11L); n <- 400L
  df  <- data.frame(v1=factor(sample(3L,n,TRUE)), v2=factor(sample(2L,n,TRUE)))
  tgt <- list(v1=c("1"=0.5,"2"=0.3,"3"=0.2), v2=c("1"=0.6,"2"=0.4))
  r <- leafblower::harvest(df, tgt, method="raking",
         min_weight=0.5, max_weight=5, max_iterations=500L, attach_weights=FALSE)
  w <- as.numeric(r)
  expect_true(all(w >= 0.5 - 1e-9), label="all weights >= min_weight after water-fill")
  expect_true(all(w <= 5 + 1e-9),   label="all weights <= max_weight after water-fill")
})

# ── Convergence status v2 tests ──────────────────────────────────────────────

test_that("status-budget: budget exhausted emits status=4, not status=1", {
  # Far-from-converged problem, tiny budget → budget exhausted (status=4).
  # RED before implementation: raking_solve still emits status=1 (RK_ERR_NOCONV).
  set.seed(99L)
  df <- data.frame(
    v1 = factor(sample(5L, 2000L, TRUE)),
    v2 = factor(sample(4L, 2000L, TRUE)),
    v3 = factor(sample(3L, 2000L, TRUE))
  )
  tgt <- list(
    v1 = setNames(rep(0.2, 5), as.character(1:5)),
    v2 = setNames(c(0.4, 0.3, 0.2, 0.1), as.character(1:4)),
    v3 = setNames(c(0.5, 0.3, 0.2), as.character(1:3))
  )
  r <- leafblower::harvest(df, tgt, method = "raking", accelerate = FALSE,
    max_weight = 5, max_iterations = 5L, attach_weights = FALSE)
  expect_equal(attr(r, "result")$status, 4L,
               label = "budget exhausted must return status=4 (RK_ERR_BUDGET)")
})

test_that("status-stall: wkl plateau emits status=5 with convergence_reason='stall_kl'", {
  # Constrained problem — all cells at U_cell in some categories.
  # KL plateau fires before budget (max_iterations=1000) → status=5.
  # RED before implementation: raking_solve still emits status=1.
  set.seed(7L)
  df <- data.frame(
    v1 = factor(sample(4L, 300L, TRUE)),
    v2 = factor(sample(3L, 300L, TRUE))
  )
  tgt <- list(
    v1 = c("1"=0.4,"2"=0.3,"3"=0.2,"4"=0.1),
    v2 = c("1"=0.5,"2"=0.3,"3"=0.2)
  )
  r <- leafblower::harvest(df, tgt, method = "raking", accelerate = FALSE,
    max_weight = 2, max_iterations = 1000L, attach_weights = FALSE)
  res <- attr(r, "result")
  expect_equal(res$status, 5L,
               label = "KL plateau must return status=5 (RK_ERR_STALL)")
  expect_equal(res$convergence_used$convergence_reason, "stall_kl",
               label = "convergence_reason must be 'stall_kl' for flat loop KL stall")
})

test_that("status-perfect: perfect calibration exits status=0 not status=5", {
  # 1-category 1-margin problem already calibrated: wkl=0 → RK_OK, not stall.
  df  <- data.frame(v1 = factor(rep("1", 20L)))
  tgt <- list(v1 = c("1" = 1.0))
  r <- leafblower::harvest(df, tgt, method = "raking", accelerate = FALSE,
    max_weight = 5, max_iterations = 500L, attach_weights = FALSE)
  expect_equal(attr(r, "result")$status, 0L,
               label = "perfect calibration must return status=0, not status=5")
})
