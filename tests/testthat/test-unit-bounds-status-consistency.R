# kxna.20 / CR-C7c: in bounds_mode="unit" the per-obs water-fill (finalize_weights_buf)
# runs AFTER the solver sets status on the pre-finalize cell iterate. A fully-pinned cell
# (feasible mass > M_cell*max_weight, forced here by heavily-skewed within-cell design
# weights) exits the water-fill clamp-only, so the RETURNED weights miss margins while the
# cell-level solve reported RK_OK. regate_unit_status() demotes that false RK_OK to STALL(5)
# when the returned proportion max_error exceeds st.tol_abs. The bug manifests in oris (its
# marginal_kl convergence rule declares RK_OK on a near-pinned cell); logit (structural
# INFEAS pre-check) and greenkhorn (greedy convergence) reach non-OK by their own guards on
# such fixtures, so the re-gate is the active demote for oris + defense-in-depth for the
# other two. The invariant "returned max_error > tol_abs => status != RK_OK" is asserted for
# all three.

# A unit-mode cell that cannot honour its target within per-obs [min,max] once the skewed
# design weights are enforced: only 15% of obs are in level "B" but its target is 0.45, and
# their design weights are heavily right-skewed, so the water-fill pins them and drifts.
.drift_fixture <- function() {
  set.seed(101L); n <- 2000L
  x <- factor(sample(c("A", "B"), n, TRUE, prob = c(0.85, 0.15)))
  y <- factor(sample(c("P", "Q"), n, TRUE, prob = c(0.5, 0.5)))
  d <- rep(1, n); d[x == "B"] <- rlnorm(sum(x == "B"), meanlog = 0, sdlog = 1.6)
  list(df = data.frame(x = x, y = y), d = d,
       tgt = list(x = c(A = 0.55, B = 0.45), y = c(P = 0.5, Q = 0.5)))
}

# A genuinely-feasible unit-mode problem (uniform design weights, loose bounds).
.feasible_fixture <- function() {
  set.seed(3L); n <- 4000L
  df <- data.frame(a = factor(sample(letters[1:3], n, TRUE)),
                   b = factor(sample(LETTERS[1:3], n, TRUE)),
                   c = factor(sample(c("M", "F"), n, TRUE)))
  tgt <- list(a = setNames(c(0.4, 0.35, 0.25), letters[1:3]),
              b = setNames(c(0.5, 0.3, 0.2), LETTERS[1:3]),
              c = c(M = 0.48, F = 0.52))
  list(df = df, tgt = tgt)
}

test_that("oris unit-mode false RK_OK on a fully-pinned cell is demoted to STALL (kxna.20)", {
  f <- .drift_fixture()
  w <- suppressWarnings(harvest(f$df, f$tgt, method = "oris", bounds_mode = "unit",
                                design_weights = f$d, max_weight = 1.5, min_weight = 0,
                                max_iterations = 500L, attach_weights = FALSE))
  r <- attr(w, "result")
  # Pre-fix: status == 0L with max_error ~ 0.156 (false RK_OK). Post-fix: demoted.
  expect_gt(r$max_error, 1e-6)
  expect_false(r$status == 0L,
               label = sprintf("oris unit fully-pinned: status=%d must not be RK_OK when max_error=%.3e",
                               r$status, r$max_error))
  expect_equal(r$status, 5L, label = "demoted to STALL (best-effort valid weights)")
})

test_that("no unit-mode solver returns RK_OK while the returned solution misses tol (kxna.20)", {
  f <- .drift_fixture()
  for (meth in c("oris", "greenkhorn", "logit")) {
    res <- tryCatch({
      w <- suppressWarnings(harvest(f$df, f$tgt, method = meth, bounds_mode = "unit",
                                    design_weights = f$d, max_weight = 1.5, min_weight = 0,
                                    max_iterations = 500L, attach_weights = FALSE))
      attr(w, "result")
    }, error = function(e) NULL)  # INFEAS harvest stop() also satisfies "no false RK_OK"
    if (!is.null(res) && is.finite(res$max_error) && res$max_error > 1e-6) {
      expect_false(res$status == 0L,
                   label = sprintf("%s unit: max_error=%.3e > tol but status=RK_OK (false OK)",
                                   meth, res$max_error))
    }
  }
})

test_that("genuinely-feasible unit-mode runs still return RK_OK (no false demote, kxna.20)", {
  f <- .feasible_fixture()
  for (meth in c("oris", "greenkhorn", "logit")) {
    w <- suppressWarnings(harvest(f$df, f$tgt, method = meth, bounds_mode = "unit",
                                  max_weight = 3, min_weight = 0.2,
                                  max_iterations = 500L, attach_weights = FALSE))
    r <- attr(w, "result")
    expect_equal(r$status, 0L,
                 label = sprintf("%s feasible unit: status=%d (max_error=%.3e) must stay RK_OK",
                                 meth, r$status, r$max_error))
    expect_lt(r$max_error, 1e-6)
  }
})
