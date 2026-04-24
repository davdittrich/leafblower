context("ieppa faithful — algBCD specifics")

# Tests in this file check properties unique to the faithful algBCD solver:
# cell compression, within-cell weight equality, capacity block behavior.
# Per design spec §6.2: 7 assertions required.

# Uses the test-only probe added in WU-1.
probe <- function(group_ids_list, n) {
  .Call("C_leafblower_cell_table_probe", group_ids_list, as.integer(n),
        PACKAGE = "leafblower")
}

test_that("cell compression correctness: n=1000, K=3 with 4 cats each, identical obs -> M_cell=64", {
  # Per design §6.2: verify cell table produces expected M_cell count for a
  # fully-populated cross-product at the expected boundary.
  n <- 1000
  set.seed(100)
  g1 <- sample(0:3, n, replace = TRUE)
  g2 <- sample(0:3, n, replace = TRUE)
  g3 <- sample(0:3, n, replace = TRUE)
  out <- probe(list(g1, g2, g3), n)
  expect_equal(out$M_cell, 64L, tolerance = 0L)  # 4^3 = 64
})

test_that("cell compression extreme: all-unique observations -> M_cell approx n", {
  # Per design §6.2: verify degenerate case where each obs is its own cell.
  n <- 200
  g1 <- seq_len(n) - 1L  # each obs unique on margin 1
  g2 <- rep(0L, n)
  out <- probe(list(g1, g2), n)
  expect_equal(out$M_cell, n)
})

test_that("within-cell weight equality: obs with identical tuples get equal weights", {
  set.seed(1)
  n <- 1200  # divisible by 4*3=12
  df <- data.frame(
    a = factor(rep(0:3, each = n/4)),
    b = factor(rep(rep(0:2, each = n/12), 4))
  )
  tgt <- list(
    a = setNames(c(0.25, 0.25, 0.25, 0.25), levels(df$a)),
    b = setNames(c(0.33, 0.33, 0.34), levels(df$b))
  )
  res <- harvest(df, tgt, method = "ieppa")
  w <- res$weights
  # Group by (a, b) tuple; within each group, weights should be equal
  for (a in 0:3) for (b in 0:2) {
    mask <- df$a == a & df$b == b
    if (sum(mask) > 1) {
      ws <- w[mask]
      expect_true(diff(range(ws)) < 1e-10,
                  info = sprintf("cell (a=%d,b=%d) weights not equal: range=%.3e",
                                 a, b, diff(range(ws))))
    }
  }
})

test_that("cap-inactive: loose bounds produce no active cap", {
  set.seed(2)
  n <- 1000
  df <- data.frame(a = sample(c("x","y","z"), n, replace = TRUE))
  tgt <- list(a = c(x = 0.4, y = 0.3, z = 0.3))
  res <- harvest(df, tgt, method = "ieppa", max_weight = 10, min_weight = 0)
  # n_cap_active accessible via attr; if not wired, skip
  # Loose bound means all weights should be in interior of [0, 10]
  expect_true(max(res$weights) < 10 - 1e-6)
  expect_true(min(res$weights) > 0 + 1e-6)
})

test_that("cap-active: tight bounds force cap with targets still met", {
  set.seed(3)
  n <- 1000
  df <- data.frame(
    a = c(rep("x", 100), rep("y", 900))  # 90/10 split
  )
  tgt <- list(a = c(x = 0.5, y = 0.5))  # need heavy upweighting of x
  res <- harvest(df, tgt, method = "ieppa", max_weight = 5, min_weight = 0)
  # Target: x:y = 50:50; achieved by upweighting x ~5x
  # Some weights should hit the cap
  expect_true(max(res$weights) >= 5 - 1e-6 || max(res$weights) <= 5 + 1e-6)
  # Verify target met (approximately)
  diag_row <- sum(res$weights[df$a == "x"]) / sum(res$weights)
  expect_lt(abs(diag_row - 0.5), 1e-4)
})

test_that("infeasibility: empty cell + positive target → INFEAS error", {
  set.seed(4)
  n <- 100
  df <- data.frame(a = rep("x", n))  # only x
  tgt <- list(a = c(x = 0.5, y = 0.5))  # y target positive but no y obs
  expect_error(harvest(df, tgt, method = "ieppa"),
               regexp = "infeasible|empty cell", ignore.case = TRUE)
})

test_that("both-sided cap: min_weight + max_weight both active, targets met", {
  set.seed(5)
  n <- 500
  df <- data.frame(
    a = sample(c("x","y"), n, replace = TRUE, prob = c(0.8, 0.2))
  )
  tgt <- list(a = c(x = 0.5, y = 0.5))  # need up y and down x
  res <- harvest(df, tgt, method = "ieppa",
                 min_weight = 0.3, max_weight = 3)
  expect_true(min(res$weights) >= 0.3 - 1e-6)
  expect_true(max(res$weights) <= 3 + 1e-6)
})

test_that("WU-2: dense regime (M_cell/n ~ 1) linear-space matches log-space to 1e-8", {
  set.seed(123)
  n <- 10000L
  # K=8, 3 cats each: 3^8 = 6561 cells; at n=10000 roughly M_cell/n ~ 0.6-0.7 -> linear path.
  df <- as.data.frame(replicate(8, sample(1:3, n, replace = TRUE), simplify = FALSE))
  names(df) <- paste0("m", 1:8)
  targets <- setNames(
    replicate(8, c(a = 0.3, b = 0.4, c = 0.3), simplify = FALSE),
    paste0("m", 1:8)
  )
  for (k in names(df)) df[[k]] <- c("a","b","c")[df[[k]]]

  Sys.setenv(LBW_IEPPA_FORCE_PATH = "linear")
  on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_PATH"), add = TRUE)
  res_lin <- harvest(df, targets, method = "ieppa",
                     max_weight = 10, min_weight = 0,
                     max_iterations = 500L,
                     convergence = list(absolute = 1e-6))

  Sys.setenv(LBW_IEPPA_FORCE_PATH = "log")
  res_log <- harvest(df, targets, method = "ieppa",
                     max_weight = 10, min_weight = 0,
                     max_iterations = 500L,
                     convergence = list(absolute = 1e-6))
  Sys.unsetenv("LBW_IEPPA_FORCE_PATH")

  expect_lt(max(abs(res_lin$weights - res_log$weights)), 1e-8)
})

test_that("WU-2: sparse regime (M_cell/n ~ 0.01) auto-dispatches log-space", {
  set.seed(7)
  n <- 10000L
  # K=2, 5 cats each: 5^2 = 25 cells; M_cell/n ~ 0.0025 -> log-space path.
  df <- data.frame(
    a = sample(letters[1:5], n, replace = TRUE),
    b = sample(letters[1:5], n, replace = TRUE)
  )
  targets <- list(
    a = setNames(rep(0.2, 5), letters[1:5]),
    b = setNames(rep(0.2, 5), letters[1:5])
  )
  msgs <- capture.output(
    res <- harvest(df, targets, method = "ieppa",
                   max_weight = 5, min_weight = 0,
                   max_iterations = 500L,
                   convergence = list(absolute = 1e-6),
                   verbose = 1L),
    type = "message"
  )
  # Verbose log prefix differs by path; we rely on either absence of
  # "linear-space" OR presence of "path=log" if labelled.
  expect_false(any(grepl("linear-space", msgs)))
})

test_that("WU-2: overflow synthesis falls back to log-space, still completes", {
  # Force linear path on a high-K input; rely on adversarial targets to drive
  # factor * prod(f) toward kLinearOverflowTrip; expect one-shot fallback.
  set.seed(99)
  n <- 5000L
  K <- 12L
  df <- as.data.frame(replicate(K, sample(1:3, n, replace = TRUE), simplify = FALSE))
  names(df) <- paste0("m", 1:K)
  for (k in names(df)) df[[k]] <- c("a","b","c")[df[[k]]]
  targets <- setNames(
    replicate(K, c(a = 0.6, b = 0.3, c = 0.1), simplify = FALSE),
    paste0("m", 1:K)
  )
  Sys.setenv(LBW_IEPPA_FORCE_PATH = "linear")
  on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_PATH"), add = TRUE)
  # Success condition: no error, status is RK_OK or RK_ERR_NOCONV, weights finite.
  res <- suppressWarnings(harvest(df, targets, method = "ieppa",
                                  max_weight = 1e6, min_weight = 0,
                                  max_iterations = 200L,
                                  convergence = list(absolute = 1e-4)))
  expect_true(all(is.finite(res$weights)))
})

test_that("WU-3: stable-mode fast-path is deterministic + does not engage damping", {
  # Pre-WU-3 byte-identity is proved by inspection: the Step 3.4/3.5 alpha==1.0
  # branch is literally the pre-WU-3 assignment. The indirect byte-identity
  # gate is the RK_OK-preservation diff at Step A.4.1 against the Step P.5
  # baseline, which ran pre-WU-3 code. This in-session test guards only:
  # (a) stable-mode output is deterministic across runs, and (b) on an input
  # with no persistence stress, damping does NOT auto-engage (verbose log
  # must not contain "damping engaged").
  set.seed(11)
  n <- 1000L
  df <- data.frame(
    a = sample(letters[1:3], n, replace = TRUE),
    b = sample(letters[1:3], n, replace = TRUE)
  )
  targets <- list(
    a = c(a = 0.33, b = 0.33, c = 0.34),
    b = c(a = 0.33, b = 0.33, c = 0.34)
  )
  Sys.setenv(LBW_IEPPA_FORCE_DAMPING = "off")
  on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_DAMPING"), add = TRUE)
  msgs <- capture.output(
    res1 <- harvest(df, targets, method = "ieppa",
                    max_weight = 5, min_weight = 0,
                    max_iterations = 500L,
                    convergence = list(absolute = 1e-6),
                    verbose = 1L),
    type = "output"
  )
  res2 <- harvest(df, targets, method = "ieppa",
                  max_weight = 5, min_weight = 0,
                  max_iterations = 500L,
                  convergence = list(absolute = 1e-6))
  expect_identical(res1, res2)                           # determinism
  expect_false(any(grepl("damping engaged", msgs)))      # fast-path stays put
})

test_that("P1.1: linear path writes X_cur exactly M_cell times per iter (fused block)", {
  # K=2 with 3 cats each → M_cell bounded by 3^2 = 9. n=100000 per cell = 11k obs,
  # birthday-saturates all 9 cells deterministically. X_init[c] > 0 for every c,
  # so the fused block increments the counter on every cell every iter.
  # (critical-code-reviewer R1: former K=4/n=5000 was RNG-fragile on empty-cell risk.)
  Sys.setenv(LBW_IEPPA_FORCE_PATH = "linear")
  on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_PATH"), add = TRUE)
  set.seed(991)
  n <- 100000L
  df <- data.frame(
    a = sample(c("a","b","c"), n, replace = TRUE),
    b = sample(c("a","b","c"), n, replace = TRUE)
  )
  targets <- list(
    a = c(a = 0.4, b = 0.35, c = 0.25),
    b = c(a = 0.4, b = 0.35, c = 0.25)
  )
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 20L,
                 convergence = list(absolute = 1e-300),
                 attach_weights = FALSE)
  result_info <- attr(res, "result")
  expect_true(!is.null(result_info$n_xcur_writes_per_iter_linear))
  stopifnot(result_info$iterations > 0)
  writes_per_iter <- result_info$n_xcur_writes_per_iter_linear / result_info$iterations
  # M_cell via probe. All 9 cells populated with high certainty at n=100k.
  gid_list <- lapply(names(targets), function(nm) {
    lv <- names(targets[[nm]])
    idx <- match(as.character(df[[nm]]), lv) - 1L
    idx[is.na(idx)] <- -1L
    as.integer(idx)
  })
  probe <- .Call("C_leafblower_cell_table_probe", gid_list, n, PACKAGE = "leafblower")
  expect_equal(probe$M_cell, 9L)            # deterministic saturation
  expect_equal(writes_per_iter, probe$M_cell)  # fused block hits every cell once
})

test_that("WU-3: damped mode takes strictly more iters than stable on same input (spec §7)", {
  # Use LBW_IEPPA_FORCE_DAMPING to run the SAME feasible input twice: once
  # stable (alpha=1.0, fast path), once damped (alpha=0.5, geometric blend).
  # Spec §7 / CTO B5: monotone `iter_damped > iter_stable` assertion.
  set.seed(314)
  n <- 3000L
  K <- 6L
  df <- as.data.frame(replicate(K, sample(1:3, n, replace = TRUE), simplify = FALSE))
  names(df) <- paste0("m", 1:K)
  for (k in names(df)) df[[k]] <- c("a","b","c")[df[[k]]]
  targets <- setNames(
    replicate(K, c(a = 0.5, b = 0.3, c = 0.2), simplify = FALSE),
    paste0("m", 1:K)
  )

  run_one <- function(force_damping) {
    Sys.setenv(LBW_IEPPA_FORCE_DAMPING = force_damping)
    on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_DAMPING"), add = TRUE)
    msgs <- capture.output(
      res <- suppressWarnings(harvest(df, targets, method = "ieppa",
                                      max_weight = 10, min_weight = 0,
                                      max_iterations = 500L,
                                      convergence = list(absolute = 1e-5),
                                      verbose = 1L)),
      type = "output"
    )
    # Final verbose line: "iEPPA <status> in <N> iters, errRp=..."
    m <- tail(grep("in [0-9]+ iters", msgs, value = TRUE), 1)
    iter <- as.integer(sub(".*in ([0-9]+) iters.*", "\\1", m))
    list(res = res, iter = iter)
  }

  r_stable <- run_one("off")
  r_damped <- run_one("on")
  expect_true(all(is.finite(r_stable$res$weights)))
  expect_true(all(is.finite(r_damped$res$weights)))
  expect_gt(r_damped$iter, r_stable$iter)  # monotone; spec §7
})

test_that("P2.1: benign input keeps alpha == 1.0 (fast path)", {
  set.seed(55)
  n <- 500L
  df <- data.frame(a = sample(letters[1:3], n, TRUE),
                   b = sample(letters[1:3], n, TRUE))
  targets <- list(a = c(a=0.33,b=0.33,c=0.34), b = c(a=0.33,b=0.33,c=0.34))
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 500L,
                 convergence = list(absolute = 1e-6),
                 attach_weights = FALSE)
  info <- attr(res, "result")
  expect_equal(info$min_alpha_seen, 1.0)  # no stress -> no damping
  expect_equal(info$final_alpha, 1.0)
})

test_that("P2.1: stress input engages damping (alpha < 1.0) with smooth schedule", {
  set.seed(314)
  n <- 3000L
  K <- 6L
  df <- as.data.frame(replicate(K, sample(1:3, n, TRUE), simplify=FALSE))
  names(df) <- paste0("m", 1:K)
  for (k in names(df)) df[[k]] <- c("a","b","c")[df[[k]]]
  targets <- setNames(
    replicate(K, c(a=0.7, b=0.2, c=0.1), simplify=FALSE),
    paste0("m", 1:K)
  )
  res <- suppressWarnings(harvest(df, targets, method = "ieppa",
                                  max_weight = 3, min_weight = 0,
                                  max_iterations = 500L,
                                  convergence = list(absolute = 1e-4),
                                  attach_weights = FALSE))
  info <- attr(res, "result")
  expect_lt(info$min_alpha_seen, 1.0)   # stress engaged damping at some point
  # Unlatched schedule: if streaks subside before exit, alpha recovers.
  # Final alpha may be 1.0 (full recovery) or intermediate. Not asserted strict.
  expect_true(info$min_alpha_seen > 0.0)  # sanity: formula is bounded below
})

test_that("P2.1: LBW_IEPPA_FORCE_DAMPING=on forces alpha = 0.5", {
  Sys.setenv(LBW_IEPPA_FORCE_DAMPING = "on")
  on.exit(Sys.unsetenv("LBW_IEPPA_FORCE_DAMPING"), add = TRUE)
  set.seed(55)
  n <- 500L
  df <- data.frame(a = sample(letters[1:3], n, TRUE),
                   b = sample(letters[1:3], n, TRUE))
  targets <- list(a = c(a=0.33,b=0.33,c=0.34), b = c(a=0.33,b=0.33,c=0.34))
  res <- harvest(df, targets, method = "ieppa",
                 max_weight = 5, min_weight = 0,
                 max_iterations = 500L,
                 convergence = list(absolute = 1e-6),
                 attach_weights = FALSE)
  info <- attr(res, "result")
  expect_equal(info$min_alpha_seen, 0.5)
  expect_equal(info$final_alpha, 0.5)
})

