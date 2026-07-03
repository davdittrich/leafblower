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

# eb79.4: TYPEOF guards on group_ids/method/start_weights. Pre-fix, a wrong
# TYPEOF here (e.g. LENGTH()/VECTOR_ELT() on a non-VECSXP, or REAL() on an
# INTSXP) can segfault the R process rather than raise a graceful error, so
# the check runs in a callr subprocess: a crash there surfaces as a callr
# error in this process instead of aborting the whole testthat suite.
run_wrong_typeof <- function(args) {
  library(leafblower)
  do.call(.Call, c(args, list(PACKAGE = "leafblower")))
}

expect_error_isolated <- function(args, regexp) {
  if (requireNamespace("callr", quietly = TRUE)) {
    expect_error(
      callr::r(run_wrong_typeof, args = list(args = args)),
      regexp
    )
  } else {
    expect_error(run_wrong_typeof(args), regexp)
  }
}

test_that("CXX.1: group_ids as non-list (atomic vector) errors gracefully (no crash)", {
  args <- make_call_args()
  args[[2]] <- as.integer(rep(0L, 12))  # wrong TYPEOF: INTSXP vector, not VECSXP list
  expect_error_isolated(args, "group_ids must be a list")
})

test_that("CXX.1: method as non-character (integer) errors gracefully (no crash)", {
  args <- make_call_args()
  args[[8]] <- 1L  # wrong TYPEOF: INTSXP, not STRSXP
  expect_error_isolated(args, "method must be a length-1 character string")
})

test_that("CXX.1: start_weights as non-numeric (integer vector) errors gracefully (no crash)", {
  args <- make_call_args()
  args[[11]] <- as.integer(rep(1L, 12))  # wrong TYPEOF: INTSXP, not REALSXP
  expect_error_isolated(args, "start_weights must be numeric")
})

# CR-D10 (j7x8.10): direct .Call bypasses R's match.arg(). An unrecognized
# method string previously validated as RAKING but executed ORIS (silent
# contract mismatch). The bridge must now reject it, naming the bad string.
test_that("CR-D10: unrecognized method string errors gracefully naming it", {
  args <- make_call_args()
  args[[8]] <- "grg"  # slot 7: method, not in kAlgMap
  expect_error(do_call_bridge(args), "unrecognized algorithm.*grg")
})

test_that("CR-D10: recognized method strings still calibrate (regression)", {
  for (m in c("auto", "oris", "oris_soft", "raking", "sinkhorn", "greg",
              "greenkhorn", "logit", "newton_kl", "chebyshev")) {
    args <- make_call_args()
    args[[8]] <- m
    expect_no_error(do_call_bridge(args))
  }
})

# CR-D9 (j7x8.9): bad scalar args (wrong type / length) now route through the
# deferred-throw path (scalar_real/scalar_int set pre_error, the single throw ->
# catch -> Rf_error site fires AFTER the RAII scope unwinds) instead of calling
# Rf_error directly mid-parse (which longjmp-leaked live std::vectors). Verify
# the graceful error still fires and names the offending parameter.
test_that("CR-D9: bad scalar arg errors gracefully via deferred throw", {
  args <- make_call_args()
  args[[6]] <- "x"  # slot 5: min_weight as character (not REALSXP)
  expect_error(do_call_bridge(args), "min_weight.*must be a length-1 numeric")

  args <- make_call_args()
  args[[7]] <- c(1.0, 2.0)  # slot 6: max_weight length 2 (not length-1)
  expect_error(do_call_bridge(args), "max_weight.*must be a length-1 numeric")

  args <- make_call_args()
  args[[9]] <- as.double(50)  # slot 8: verbose (scalar_int) given a double
  expect_error(do_call_bridge(args), "verbose.*must be a length-1 integer")
})
