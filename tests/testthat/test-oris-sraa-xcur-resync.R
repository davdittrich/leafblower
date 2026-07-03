# ──────────────────────────────────────────────────────────────────────────────
# CR-A4 (mxcl.4): ORIS SRAA convergence / best-iterate must evaluate the ACCEPTED
# iterate, not a stale pre-sweep X_cur.
#
# On the LOG path (use_linear = n/M_cell < 2 is FALSE, i.e. n/M_cell >= 2), the
# hoisted f_eval_lf refreshes X_tilde but leaves X_cur == unpack of the PRE-sweep
# input (oris.cpp:762 / 794-796). The SRAA convergence + best-iterate check reads
# X_cur, so it evaluates the iterate one BCD sweep behind the accepted lf_flat.
# (Separately, on the LINEAR path X_cur is maintained in place and is correct
# EXCEPT after a rejected AA step — sraa.hpp:226 leaves the globals at X_AA while
# lf_flat reverts to the plain step.)
#
# Observable symptom: res$metric_first_check (the reported first-check metric) is
# pinned at the stale X_cur == X_init value instead of the accepted first iterate.
# On a log-path fixture of two INDEPENDENT margins, one BCD sweep converges each
# margin near-exactly, so the accepted first iterate is ~machine-zero while X_init
# (uniform weights) is far from target. Measured pre-fix: metric_first_check
# equals X_init's marginal_kl to 5 digits (7.5e-2), while the converged weights
# have marginal_kl ~1e-16.
# ──────────────────────────────────────────────────────────────────────────────

test_that("ORIS SRAA log-path check evaluates accepted iterate, not stale X_cur (CR-A4)", {
  set.seed(59)
  n  <- 300L
  # 3x3 = up to 9 cells at n=300 => n/M_cell ~ 33 >= 2 => LOG path.
  df <- data.frame(
    a = factor(sample(letters[1:3], n, replace = TRUE)),
    b = factor(sample(LETTERS[1:3], n, replace = TRUE))
  )
  mk <- function(nm) { v <- runif(3, 0.5, 1.5); v <- v / sum(v); setNames(v, nm) }
  target <- list(a = mk(letters[1:3]), b = mk(LETTERS[1:3]))

  # Fresh, package-independent max-over-margins marginal_kl (matches the SRAA
  # convergence metric = select_metric(MARGINAL_KL→KL) = CellMetrics.kl).
  marg_kl <- function(w) {
    W <- sum(w); eps <- 1e-10; kl <- 0
    for (v in names(target)) {
      Sp <- tapply(w, df[[v]], sum)[names(target[[v]])] / W
      T  <- target[[v]]
      kl <- max(kl, sum(ifelse(T > 0, T * log((T + eps) / (Sp + eps)), 0)))
    }
    kl
  }
  # oris X_init proportions == uniform-weight (design) proportions.
  X_init_kl <- marg_kl(rep(1, n))

  w <- suppressWarnings(harvest(df, target, method = "oris", accelerate = TRUE,
                                max_iterations = 2000L, attach_weights = FALSE))
  r <- attr(w, "result")
  expect_identical(r$status, 0L)

  # Primary: the reported first-check metric must reflect the accepted first
  # iterate (a real converging value, here ~1.6e-3), NOT the stale pre-sweep
  # X_cur == X_init. Pre-fix this equals X_init_kl (7.5e-2) exactly and fails;
  # post-fix it drops an order of magnitude below X_init.
  expect_lt(r$metric_first_check, X_init_kl / 10)

  # DoD invariant on the LOG path: the reported max_error matches the max margin
  # residual independently recomputed from the RETURNED weights (both realize the
  # accepted, clamped iterate). This is the CR-A5-class gate exercised on the SRAA
  # log path specifically.
  max_err_ret <- {
    W <- sum(w); e <- 0
    for (v in names(target))
      e <- max(e, max(abs(tapply(w, df[[v]], sum)[names(target[[v]])] / W - target[[v]])))
    e
  }
  expect_equal(max_err_ret, r$max_error, tolerance = 1e-8)
})
