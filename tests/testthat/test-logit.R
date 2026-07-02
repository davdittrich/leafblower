test_that("F(0) = 1 for logit link", {
  expect_equal(.Call("C_logit_F_at_zero", 0.5, 5.0), 1.0, tolerance = 1e-12)
})

test_that("F(u) stays in [L,U] for logit link", {
  vals <- .Call("C_logit_range_check", 0.5, 5.0, as.double(seq(-10, 10, by=0.5)))
  expect_true(all(vals >= 0.5 - 1e-12))
  expect_true(all(vals <= 5.0 + 1e-12))
})

test_that("H prime equals F for logit link (numerical diff)", {
  result <- .Call("C_logit_Hprime_check", 0.5, 5.0, 1.0)
  expect_equal(result, 0.0, tolerance = 1e-8)
})

test_that("exp link: F(u) = exp(u)", {
  expect_equal(.Call("C_logit_F_at_zero", 0.0, Inf), 1.0, tolerance = 1e-12)
})

test_that("H'(u) = F(u) holds near safe_exp clamp boundary (u=559 for L=0,U=5)", {
  # logit_scale = (5-0)/((5-1)*(1-0)) = 1.25
  # logit_scale * 559 = 698.75 (just below the safe_exp clamp at 700)
  # safe_exp clamp preserves H'(u) = F(u) because both F and H use the same
  # safe_exp(logit_scale*u) value; the algebraic identity is maintained.
  result <- .Call("C_logit_Hprime_check", 0.0, 5.0, 559.0)
  expect_equal(result, 0.0, tolerance = 1e-4)  # wider tol: finite-diff truncation error near saturation (F(u)≈U=5)
})

# ---------------------------------------------------------------------------
# eb79.3 (T3): signed-Δz Newton step-norm guard (logit_calib.cpp:238-249)
#
# The step-size cap now accumulates the SIGNED per-cell logit shift
# Δz_c = Σ_k Δλ[k,g] and takes |·| once, matching the actual (signed) cell
# coordinate z_c = Σ_k λ[k,g]. The prior abs-per-margin sum over-estimated
# |Δz_c| (triangle inequality) and throttled alpha on every K>=2 cell.
#
# INVARIANT (build-independent): the fix only changes throttling, NOT the
# fixed point. A converged logit solution must (a) recover every margin
# target, (b) have finite weights inside [min_weight, max_weight]. On the
# collinear K>=6 fixture, (b) is also the observable guard for the agy RISK:
# dual-λ null-space drift / catastrophic cancellation in z_c would surface as
# non-finite weights or bound violations. The dual λ itself is not exposed to
# R, so max|λ| non-inflation is verified in C++ by the orchestrator; the R
# layer asserts the visible symptom is absent.
#
# CROSS-VERSION assertions (identical weights vs pre-fix, iterations not
# worse, max|λ| not inflated) require running BOTH pre-fix and post-fix
# builds; those reference numbers are filled by the orchestrator (constants
# flagged with ORCH_FILL below). Do NOT self-certify pass/fail on them.
# ---------------------------------------------------------------------------

test_that("eb79.3: K>=2 logit converges to a valid fixed point (targets recovered, in-bounds)", {
  library(leafblower)
  set.seed(11)
  n <- 300
  data <- data.frame(
    a = factor(sample(c("1", "2", "3"), n, replace = TRUE)),
    b = factor(sample(c("x", "y"), n, replace = TRUE))
  )
  target <- list(a = c("1" = 1/3, "2" = 1/3, "3" = 1/3), b = c(x = 0.5, y = 0.5))
  w <- harvest(data, target, method = "logit", min_weight = 0.2, max_weight = 5,
               max_iterations = 500L, attach_weights = FALSE)
  r <- attr(w, "result")

  expect_equal(r$status, 0L, label = "status must be RK_OK=0")
  expect_true(all(is.finite(w)), label = "all weights finite (no λ blowup)")
  expect_true(all(w >= 0.2 - 1e-9 & w <= 5 + 1e-9), label = "weights within bounds")

  # Fixed-point invariant: weighted margins recover targets * n.
  expect_equal(sum(w[data$a == "1"]), n/3, tolerance = 1e-4)
  expect_equal(sum(w[data$b == "x"]), n/2, tolerance = 1e-4)
  expect_equal(sum(w), n, tolerance = 1e-6)

  # Regression guard: signed-Δz must not increase iterations vs the pre-fix
  # abs-sum guard. Orchestrator measured both builds on this fixture (seed 11):
  # pre-fix iters=14, post-fix iters=14 (iteration-neutral here — abs vs signed
  # only differ on cells with >=2 opposing-sign active margins; this 2-margin
  # fixture rarely hits that). Hard constant so any future regression trips it.
  ITERS_PREFIX_K2 <- 14L
  expect_lte(r$iterations, ITERS_PREFIX_K2)
})

test_that("eb79.3: collinear K>=6 fixture — no λ null-space blowup (agy RISK)", {
  library(leafblower)
  set.seed(23)
  n <- 200
  base <- factor(c(rep("A", n/2), rep("B", n/2)))
  base2 <- factor(rep(c("P", "Q", "P", "Q"), length.out = n))  # 50/50, orthogonal-ish to base
  # 6 margins, deliberately rank-deficient: v1..v3 are identical copies of base
  # (perfectly collinear), v4..v6 identical copies of base2 → AA^T is rank-2 but
  # embedded in a 12-column dual space, i.e. a large null space. This is exactly
  # the setting the signed-Δz cap could expose to dual drift.
  data <- data.frame(
    v1 = base, v2 = base, v3 = base,
    v4 = base2, v5 = base2, v6 = base2
  )
  tgt_ab <- c(A = 0.5, B = 0.5)
  tgt_pq <- c(P = 0.5, Q = 0.5)
  target <- list(v1 = tgt_ab, v2 = tgt_ab, v3 = tgt_ab,
                 v4 = tgt_pq, v5 = tgt_pq, v6 = tgt_pq)

  w <- harvest(data, target, method = "logit", min_weight = 0.1, max_weight = 10,
               max_iterations = 500L, attach_weights = FALSE)
  r <- attr(w, "result")

  # (b) Observable guard for the agy RISK: finite, in-bounds weights.
  expect_true(all(is.finite(w)),
              label = "collinear fixture: weights finite (no catastrophic cancellation in z_c)")
  expect_true(all(w >= 0.1 - 1e-9 & w <= 10 + 1e-9),
              label = "collinear fixture: weights within bounds (no null-space λ drift blowup)")

  # This fully-collinear fixture is DEGENERATE for the logit link: the solver
  # settles to all weights at min_weight (z_c -> -inf) in 1 iteration, which
  # preserves each margin as a PROPORTION (0.5/0.5) but leaves Sum(w) != n. That
  # is a pre-existing logit behavior (tracked separately — see epic eb79), NOT a
  # product of the signed-Δz change, so we do NOT assert Sum(w)==n here.
  #
  # eb79.16: collinear result was always infeasible; status now honestly reports
  # STALL. Pre-fix the scale-blind proportion metric certified this
  # rank-deficient stall as status=0/max_error=0. The absolute-count convergence
  # gate now refuses to declare success on a violated absolute constraint (the
  # true margin residual is ~90 counts, ~0.45 in proportion units), so status is
  # honestly STALL/INFEAS (!= 0). The finite + in-bounds guards above still hold.
  expect_true(r$status != 0L,
    label = "collinear fixture: honest STALL/INFEAS (was misreported status=0 pre-eb79.16)")
  # The honest reported error reflects the true absolute-count violation, not 0.
  expect_true(r$max_error > 0.1,
    label = "collinear fixture: reported max_error reflects real violation (~0.45), not 0")
})

# ---------------------------------------------------------------------------
# eb79.16 (E1): scale-aware logit convergence + absolute-space reported error.
# The convergence gate and residual-class reported fields (max_error,
# best_error, mean_error, grake_norm) are computed against the true target
# total n, not the current weight sum. kl/chi2 stay proportion-space.
# ---------------------------------------------------------------------------

test_that("eb79.16: collinear fixture reports honest infeasibility (status != 0, max_error > 0.1)", {
  library(leafblower)
  set.seed(23)
  n <- 200
  base <- factor(c(rep("A", n/2), rep("B", n/2)))
  base2 <- factor(rep(c("P", "Q", "P", "Q"), length.out = n))
  data <- data.frame(
    v1 = base, v2 = base, v3 = base,
    v4 = base2, v5 = base2, v6 = base2
  )
  tgt_ab <- c(A = 0.5, B = 0.5)
  tgt_pq <- c(P = 0.5, Q = 0.5)
  target <- list(v1 = tgt_ab, v2 = tgt_ab, v3 = tgt_ab,
                 v4 = tgt_pq, v5 = tgt_pq, v6 = tgt_pq)

  w <- harvest(data, target, method = "logit", min_weight = 0.1, max_weight = 10,
               max_iterations = 500L, attach_weights = FALSE)
  r <- attr(w, "result")

  expect_true(r$status != 0L, label = "collinear: honest STALL/INFEAS")
  expect_true(r$max_error > 0.1,
              label = "collinear: absolute-space max_error ~0.45, not scale-blind 0")
})

test_that("eb79.16: non-default convergence cfg on a feasible fixture still converges (no spurious STALL)", {
  library(leafblower)
  set.seed(11)
  n <- 300
  data <- data.frame(
    a = factor(sample(c("1", "2", "3"), n, replace = TRUE)),
    b = factor(sample(c("x", "y"), n, replace = TRUE))
  )
  target <- list(a = c("1" = 1/3, "2" = 1/3, "3" = 1/3), b = c(x = 0.5, y = 0.5))

  # metric = "kl": absolute gate must not flip a legitimately-converged case.
  w_kl <- harvest(data, target, method = "logit", min_weight = 0.2, max_weight = 5,
                  max_iterations = 500L, attach_weights = FALSE,
                  convergence = list(metric = "kl"))
  r_kl <- attr(w_kl, "result")
  expect_equal(r_kl$status, 0L, label = "metric=kl feasible fixture: status RK_OK=0")
  expect_equal(sum(w_kl), n, tolerance = 1e-6)

  # absolute shorthand (metric="max_err", rule="threshold"): the absolute gate
  # scales with the user tol and must NOT spuriously flip a feasible case.
  w_thr <- harvest(data, target, method = "logit", min_weight = 0.2, max_weight = 5,
                   max_iterations = 500L, attach_weights = FALSE,
                   convergence = list(absolute = 1e-3))
  r_thr <- attr(w_thr, "result")
  expect_equal(r_thr$status, 0L, label = "absolute-rule feasible fixture: status RK_OK=0")
  # Loose user tol (absolute=1e-3): the absolute gate stops at max_abs_resid<=1e-3*n,
  # so Sum(w) is within the user's chosen tolerance of n (NOT machine-exact) — the
  # point of this test is that a feasible case is NOT spuriously flipped to STALL.
  expect_equal(sum(w_thr), n, tolerance = 1e-2)
})

test_that("eb79.16: tight-but-feasible fixture converges (status=0, Sum(w)=n)", {
  library(leafblower)
  set.seed(11)
  n <- 300
  data <- data.frame(
    a = factor(sample(c("1", "2", "3"), n, replace = TRUE)),
    b = factor(sample(c("x", "y"), n, replace = TRUE))
  )
  target <- list(a = c("1" = 1/3, "2" = 1/3, "3" = 1/3), b = c(x = 0.5, y = 0.5))
  w <- harvest(data, target, method = "logit", min_weight = 0.2, max_weight = 5,
               max_iterations = 500L, attach_weights = FALSE)
  r <- attr(w, "result")
  expect_equal(r$status, 0L, label = "tight-but-feasible: status RK_OK=0")
  expect_equal(sum(w), n, tolerance = 1e-6)
})

test_that("eb79.20: logit populates l1_weight_change = mean(|w - 1|) (was struct default)", {
  library(leafblower)
  # Deterministic skewed-but-feasible fixture (no RNG) — bit-identical to the Python
  # parity test (test_python.py test_logit_l1_and_objective_populated_eb79_20): data
  # 60/40, targets 45/55 → calibration moves weights off the design value of 1, so
  # l1_weight_change is strictly positive and the hand-check below is real (not 0==0).
  n <- 400
  data <- data.frame(
    a = factor(c(rep("x", 240), rep("y", 160))),
    b = factor(ifelse((seq_len(n) - 1L) %% 100L < 60L, "p", "q"))
  )
  target <- list(a = c(x = 0.45, y = 0.55), b = c(p = 0.45, q = 0.55))
  w <- harvest(data, target, method = "logit", min_weight = 0.2, max_weight = 5,
               max_iterations = 500L, attach_weights = FALSE)
  r <- attr(w, "result")
  expect_equal(r$status, 0L, label = "skewed-feasible: status RK_OK=0")

  # l1_weight_change = Σ_i|Δw|/n. Default design weights are 1 (Σw=n), so the
  # calibrated-minus-design L1 equals mean(|w - 1|). Was left at 0 (struct default)
  # before eb79.20; now populated in logit_calib.cpp exit block.
  expect_true(is.finite(r$l1_weight_change))
  expect_gt(r$l1_weight_change, 0)                                   # weights actually moved
  expect_equal(r$l1_weight_change, mean(abs(w - 1)), tolerance = 1e-6)
})
