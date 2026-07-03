# ──────────────────────────────────────────────────────────────────────────────
# CR-A5 (mxcl.5): cross-solver returned-weights ≡ reported-metrics regression gate.
#
# The three Blocking defects in epic CR-A share one root symptom: the weights
# returned to the caller do NOT realize the metrics reported in the result struct.
#   • CR-A1 (mxcl.1): Cholesky triangle mismatch → greg/logit reported converged
#     while returning a diagonal-solve iterate (1-step err 3.7e-2).
#   • CR-A2 (mxcl.2): chebyshev warm-start pairing → returned ≈design weights
#     (err 0.215) while metrics claimed calibrated.
#   • CR-A6 (mxcl.6): ORIS SRAA reported pre-clamp max_error (2e-16) while the
#     returned clamped weights were 26% off.
#
# No prior test independently re-derives metrics from the RETURNED weights across
# ALL solvers — each suite trusts its own self-reported numbers. This gate closes
# that gap permanently: for every method it recomputes max_error / kl / chi2 from
# the returned weights using a FRESH base-R implementation (never the package's own
# metric code) and asserts they match the reported fields.
#
# Formulas mirror compute_cell_metrics (src/calib_dispatch.hpp:239), kMetricEps=1e-10:
#   S_p    = Σ_{i∈cat(k,j)} w_i / W        (achieved proportion, W = Σw)
#   maxerr = max_{k,j} |S_p − T_kj|                              (errRp, MAX)
#   kl     = max_k Σ_j [T>0] T·log((T+ε)/(S_p+ε))               (reverse-KL, MAX)
#   chi2   = Σ_{k,j: pop>ε} (obs − pop)² / pop,  obs=Σw, pop=T·W (Pearson)
#
# Per-method tolerances are empirically derived (sweep 2026-07-03, n=600 fixture)
# with generous safety headroom. max_error/kl are scale-invariant proportions
# (tight); chi2 is W-scaled (looser). The epic's bugs were 0.037–0.26 — every
# tolerance below leaves ≥1e3× margin, so the gate never false-fails yet catches
# any regression of this class. Observed |Δ| in the sweep:
#   raking/sinkhorn/greg/chebyshev/oris/oris_soft/auto : ~1e-16 (FP reassociation)
#   newton_kl : Δkl 3.9e-14, Δchi2 4.7e-11 (BLAS reassociation)
#   greenkhorn: Δmaxerr 8e-9 (entropic, report-vs-final iterate drift; Σw≠n)
#   logit     : Δmaxerr 2.7e-7 (best-iterate/final-Newton drift; see CR-C7/CR-C8)
#
# NOTE: CR-A2 (mxcl.2) has LANDED (commit 5fdafb7), so chebyshev is asserted GREEN
# like every other method — the ticket's "mark chebyshev expected-RED" provision is
# obsolete now that its dependency is closed.
# ──────────────────────────────────────────────────────────────────────────────

test_that("returned weights realize the reported max_error/kl/chi2 for every solver (CR-A5)", {
  set.seed(99)
  n   <- 600L
  df  <- data.frame(
    x = factor(sample(c("a", "b", "c"), n, replace = TRUE)),
    y = factor(sample(c("p", "q"),      n, replace = TRUE))
  )
  target <- list(x = c(a = 1/3, b = 1/3, c = 1/3), y = c(p = 0.5, q = 0.5))

  # Fresh, package-independent re-derivation of the three cell metrics from
  # returned weights. Uses only base-R tapply/log — no leafblower internals.
  rederive <- function(w) {
    W <- sum(w); eps <- 1e-10
    maxerr <- 0; kl <- 0; chi2 <- 0
    for (v in names(target)) {
      obs <- tapply(w, df[[v]], sum)[names(target[[v]])]
      T   <- target[[v]]
      Sp  <- obs / W
      maxerr <- max(maxerr, max(abs(Sp - T)))
      kl_k   <- sum(ifelse(T > 0, T * log((T + eps) / (Sp + eps)), 0))
      kl     <- max(kl, kl_k)
      pop    <- T * W
      chi2   <- chi2 + sum(ifelse(pop > eps, (obs - pop)^2 / pop, 0))
    }
    c(max_error = maxerr, kl = kl, chi2 = chi2)
  }

  # atol_maxerr, atol_kl, atol_chi2 per method (empirical sweep × safety headroom).
  tol <- list(
    raking     = c(1e-8, 1e-8, 1e-5),
    sinkhorn   = c(1e-8, 1e-8, 1e-5),
    greenkhorn = c(1e-5, 1e-8, 1e-3),   # entropic: report-vs-final drift, Σw≠n
    logit      = c(1e-5, 1e-8, 1e-5),   # best-iterate/final-Newton drift
    greg       = c(1e-8, 1e-8, 1e-5),
    chebyshev  = c(1e-8, 1e-8, 1e-5),   # CR-A2 landed → GREEN
    newton_kl  = c(1e-8, 1e-8, 1e-5),   # BLAS reassociation
    oris       = c(1e-8, 1e-8, 1e-5),
    oris_soft  = c(1e-8, 1e-8, 1e-5),
    auto       = c(1e-8, 1e-8, 1e-5)
  )

  for (m in names(tol)) {
    w <- suppressWarnings(
      harvest(df, target, method = m, max_iterations = 1000L, attach_weights = FALSE)
    )
    r  <- attr(w, "result")
    expect_identical(r$status, 0L, label = paste("status for", m))
    rd <- rederive(w)
    t  <- tol[[m]]
    expect_lt(abs(rd["max_error"] - r$max_error), t[1],
              label = paste("returned-weights max_error mismatch for", m))
    expect_lt(abs(rd["kl"]        - r$kl),        t[2],
              label = paste("returned-weights kl mismatch for", m))
    expect_lt(abs(rd["chi2"]      - r$chi2),      t[3],
              label = paste("returned-weights chi2 mismatch for", m))
  }
})
