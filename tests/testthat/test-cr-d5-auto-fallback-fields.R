# CR-D5 (leafblower-j7x8.5): on method="auto" fallback (ORIS primary NOCONV/BUDGET
# → newton_kl), the exported result must carry the winning (newton) solver's
# diagnostics, NOT stale ORIS-only fields left over from the abandoned primary run.

test_that("auto-fallback resets stale ORIS diagnostics to non-ORIS defaults (CR-D5)", {
  set.seed(11); n <- 500
  df <- data.frame(a  = factor(sample(c("x", "y", "z"), n, TRUE, c(.85, .10, .05))),
                   b  = factor(sample(c("p", "q"),      n, TRUE, c(.75, .25))),
                   cc = factor(sample(c("m", "o"),      n, TRUE, c(.6, .4))))
  tgt <- list(a = c(x = 0.34, y = 0.33, z = 0.33), b = c(p = 0.5, q = 0.5), cc = c(m = 0.5, o = 0.5))

  # Sanity: an ORIS+accelerate run on this fixture populates ORIS-only diagnostics
  # (aa_accepted_count > 0 and homotopy_levels_used >= 1), so the reset below is
  # load-bearing — without it, the fallback result would surface these stale.
  oris_primary <- attr(suppressWarnings(
    harvest(df, tgt, method = "oris", accelerate = TRUE, max_iter = 8)), "result")
  expect_gt(oris_primary$aa_accepted_count, 0L)
  expect_gte(oris_primary$homotopy_levels_used, 1L)

  # method="auto" with accelerate + tiny budget: the ORIS primary runs SRAA then
  # NOCONV/BUDGETs, triggering the newton_kl fallback.
  L <- attr(suppressWarnings(
    harvest(df, tgt, method = "auto", accelerate = TRUE, max_iter = 8)), "result")
  expect_equal(L$algorithm_used, "newton_kl")        # fallback fired

  # ORIS-only diagnostics reset to documented non-ORIS defaults (leafblower.h):
  expect_equal(L$aa_accepted_count, 0L)              # SRAA accept count (was >0 in primary)
  expect_equal(L$homotopy_levels_used, 0L)           # ORIS sets >=1; newton owns none
  expect_equal(L$n_xcur_writes_per_iter_last, 0L)
  expect_equal(L$min_alpha_seen, 1.0)
  expect_equal(L$final_alpha, 1.0)
  expect_equal(L$homotopy_final_factor, 1.0)
  expect_equal(L$greedy_sweeps_taken, 0L)
  expect_equal(L$eta_final, 0.0)
  expect_false(isTRUE(L$sraa_demoted))
  expect_equal(L$sor$omega_mean, 1.0)
  expect_equal(L$sor$min_omega, 1.0)
})
