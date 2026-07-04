# ──────────────────────────────────────────────────────────────────────────────
# jy0m.1 (CR-C5b.1): raking box-feasibility certificate.
#
# A per-margin BOX-FEASIBILITY CERTIFICATE at the raking feasibility gate: margin
# (k,j) is PROVABLY infeasible iff its target mass T_kj·n lies outside the
# achievable band [Σ L_cell, Σ U_cell] over its cells. A box-infeasible margin
# denies RK_OK even when errRp ≤ 1pp — the run falls through to the STALL detector
# (constrained optimum, status=5) instead of a false RK_OK. Rigorous NECESSARY
# condition with zero fitted threshold; silent (byte-identical) on feasible inputs.
# ──────────────────────────────────────────────────────────────────────────────

test_that("box-infeasible rare-category corner denies false RK_OK → STALL(5) (jy0m.1)", {
  # Rare category 'r' ≈1.4% empirical, target 0.005, min_weight=0.7. Best achievable
  # r-mass is min_weight·n_r ≈ 0.7·0.014 ≈ 0.0098 > target 0.005 ⇒ the r-margin is
  # box-infeasible. errRp for that margin (|0.0098-0.005| ≈ 0.0048) sits INSIDE the
  # 1pp bar, so the errRp-only gate granted a false RK_OK; the box certificate denies
  # it. This is the exact corner from the jy0m.1 evidence.
  set.seed(11); n <- 20000L
  df <- data.frame(
    a = factor(sample(c("p", "q", "r"), n, TRUE, prob = c(0.5, 0.486, 0.014))),
    b = factor(sample(c("s", "t"), n, TRUE)))
  tg <- list(a = c(p = 0.4, q = 0.595, r = 0.005), b = c(s = 0.5, t = 0.5))
  w <- suppressWarnings(harvest(df, tg, method = "raking",
                                min_weight = 0.7, max_weight = 1.3,
                                max_iterations = 500L, attach_weights = FALSE))
  r <- attr(w, "result")
  expect_equal(r$status, 5L)          # box-cert → STALL, not false RK_OK
})

test_that("feasible problem at default bounds is silent → RK_OK (jy0m.1 byte-identical)", {
  # Same rare category, but default bounds (min=0, max=Inf): every margin is
  # box-feasible (r can reach 0.005), so the scan is silent and RK_OK is preserved.
  set.seed(11); n <- 20000L
  df <- data.frame(
    a = factor(sample(c("p", "q", "r"), n, TRUE, prob = c(0.5, 0.486, 0.014))),
    b = factor(sample(c("s", "t"), n, TRUE)))
  tg <- list(a = c(p = 0.4, q = 0.595, r = 0.005), b = c(s = 0.5, t = 0.5))
  w <- suppressWarnings(harvest(df, tg, method = "raking",
                                max_iterations = 500L, attach_weights = FALSE))
  r <- attr(w, "result")
  expect_equal(r$status, 0L)          # feasible ⇒ RK_OK
  expect_lt(r$max_error, 1e-2)
})

test_that("tight-but-satisfiable bounds converge → RK_OK (jy0m.1 no over-rejection)", {
  # Near-uniform empirical + uniform targets under tight bounds [0.7, 1.5]: every
  # margin's target ratio ≈ 1 ∈ [0.7, 1.5] ⇒ box-feasible; must still converge RK_OK.
  set.seed(1); n <- 2000L
  df <- data.frame(a = factor(sample(c("p", "q", "r"), n, TRUE)),
                   b = factor(sample(c("s", "t"), n, TRUE)))
  tg <- list(a = c(p = 1/3, q = 1/3, r = 1/3), b = c(s = .5, t = .5))
  w <- suppressWarnings(harvest(df, tg, method = "raking",
                                min_weight = 0.7, max_weight = 1.5,
                                max_iterations = 500L, attach_weights = FALSE))
  r <- attr(w, "result")
  expect_equal(r$status, 0L)          # tight but satisfiable ⇒ RK_OK, not STALL
  expect_lt(r$max_error, 1e-2)
})
