library(leafblower)

# CXX.1 (leafblower-5fm8.1): r_bridge length checks vs K.
# K is derived from LENGTH(group_ids). A direct .Call with cat_counts /
# targets shorter than K (or a target vector shorter than its cat_count)
# previously OOB-read. The bridge must now raise a graceful R error.

# Build the exact 39-arg C_rk_calibrate payload for a small valid problem,
# then return the arg list so individual tests can corrupt one slot.
make_call_args <- function() {
  cats <- c("a", "b", "c")
  K <- 3L; n <- 12L
  # 0-indexed group ids per margin (deterministic, all levels present)
  g <- lapply(seq_len(K), function(k) as.integer(((seq_len(n) - 1L) + k) %% 3L))
  cat_counts_r <- rep(3L, K)
  targets_r <- replicate(K, rep(1 / 3, 3L), simplify = FALSE)
  list(
    "C_rk_calibrate",
    g,                       # 1: group_ids  -> K = length(g)
    cat_counts_r,            # 2: cat_counts
    targets_r,               # 3: targets
    as.integer(n),           # 4: n_obs
    as.double(0),            # 5: min_weight
    as.double(100),          # 6: max_weight
    as.character("raking"),  # 7: method
    as.integer(0L),          # 8: verbose
    as.integer(50L),         # 9: max_iterations
    as.double(rep(1, n)),    # 10: start_weights
    -1.0, -1.0,              # 11-12: capacity_penalty, alm_penalty
    as.double(1e-6),         # 13: tol_abs
    0L,                      # 14: bounds_mode
    1L,                      # 15: homotopy_levels
    1.0, 1.0, 1.0,           # 16-18: homotopy start/end/budget
    as.character("none"),    # 19: scheduler
    as.character("none"),    # 20: eta_schedule
    1.0, 1.0, 1.0,           # 21-23: eta start/end/power
    as.double(0.01),         # 24: pct_tol
    as.double(1e-6),         # 25: absolute_tol
    1L, 1L, 1L,              # 26-28: metric, rule, stop_when
    0L, 0L,                  # 29-30: sor_enabled, sor_auto
    1.5, 1.0, 1.5, 1.0, 0L, # 31-35: sor omega_init/min/max/fixed, burnin
    2L,                      # 36: sor_omega_mode_id (2=spectral)
    0L,                      # 37: accelerate
    as.double(1e-8),         # 38: newton_tsvd_ratio
    as.double(0.0)           # 39: ridge_lambda
  )
}

do_call_bridge <- function(args) {
  do.call(.Call, c(args, list(PACKAGE = "leafblower")))
}

test_that("CXX.1: baseline .Call payload calibrates without error (sanity)", {
  args <- make_call_args()
  expect_no_error(do_call_bridge(args))
})

test_that("CXX.1: cat_counts shorter than K errors gracefully (no OOB)", {
  args <- make_call_args()
  args[[3]] <- rep(3L, 2L)  # slot 2 (arg index 3): cat_counts length 2 < K=3
  expect_error(do_call_bridge(args), "cat_counts length")
})

test_that("CXX.1: targets list shorter than K errors gracefully (no OOB)", {
  args <- make_call_args()
  args[[4]] <- replicate(2L, rep(1 / 3, 3L), simplify = FALSE)  # length 2 < K=3
  expect_error(do_call_bridge(args), "targets must be a list")
})

test_that("CXX.1: target vector shorter than cat_counts[k] errors gracefully (no OOB)", {
  args <- make_call_args()
  tg <- args[[4]]
  tg[[2]] <- c(0.5, 0.5)  # length 2 != cat_counts[2] == 3
  args[[4]] <- tg
  expect_error(do_call_bridge(args), "length != cat_counts")
})
