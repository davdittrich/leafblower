# benchmarks/algo_selection_benchmark.R
# Bayesian Level Set Estimation benchmark for iEPPA vs L-BFGS-B algorithm selection.
# See: docs/superpowers/specs/2026-04-20-algo-selection-design.md

# Enforce single-threaded execution before loading any library that might
# initialise an OpenMP thread pool. OMP reads OMP_NUM_THREADS at pool creation
# time (lazy, first parallel region). Setting it here — before library() calls —
# guarantees the env var is visible at that point, eliminating CPU-contention
# confound from the timing comparisons.
Sys.setenv(OMP_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(leafblower)
  library(DiceKriging)
  library(lhs)
  library(ggplot2)
})

# Guard: set .BENCH_SOURCED <- TRUE before source()ing this file in tests
# to prevent run_benchmark() from executing.
.BENCH_SOURCED <- exists(".BENCH_SOURCED", envir = .GlobalEnv, inherits = FALSE)

# Input space bounds
BENCH_X1_RANGE <- c(4.0, 7.7)   # log10(complexity): 10K to 50M
BENCH_X2_RANGE <- c(-6.0, -3.0) # log10(tol_abs):    1e-6 to 1e-3
BENCH_THRESHOLD <- log(1.2)      # log(1.2) ≈ 0.182; L-BFGS-B wins above this

# ── bench_seed ────────────────────────────────────────────────────────────────
# Deterministic integer seed from (log_complexity, log_tol).
# Valid for log_complexity in [4, 7.7] and log_tol in [-6, -3].
# Formula produces seeds in [4e8, 7.7e8+5999], all within 32-bit integer range.
bench_seed <- function(log_complexity, log_tol) {
  a <- as.integer(round(log_complexity * 1e4))
  b <- as.integer(round(-log_tol * 1e4)) %% 10000L
  (a * 10000L) + b
}

# ── make_bench_data ───────────────────────────────────────────────────────────
# Generates n-row survey data with K categorical margins of cats_per_margin levels.
# Population proportions: Dirichlet(1,...,1) via normalised Exp(1) draws.
# Sample proportions: population * |Normal(1, 0.1)| noise, renormalised (~10% bias).
make_bench_data <- function(n, K, cats_per_margin) {
  stopifnot(n >= 1L, K >= 1L, cats_per_margin >= 2L)
  col_names <- paste0("m", seq_len(K))
  lvl_names <- lapply(seq_len(K), function(k) paste0("c", seq_len(cats_per_margin)))

  pop_props  <- lapply(seq_len(K), function(k) {
    x <- rexp(cats_per_margin); x / sum(x)
  })
  samp_props <- lapply(pop_props, function(p) {
    q <- p * abs(rnorm(length(p), mean = 1, sd = 0.1)); q / sum(q)
  })

  df_cols <- lapply(seq_len(K), function(k) {
    factor(sample(lvl_names[[k]], n, replace = TRUE, prob = samp_props[[k]]),
           levels = lvl_names[[k]])
  })
  df <- as.data.frame(setNames(df_cols, col_names))

  targets <- setNames(
    lapply(seq_len(K), function(k) setNames(pop_props[[k]], lvl_names[[k]])),
    col_names
  )
  list(df = df, targets = targets)
}

# ── time_cell ─────────────────────────────────────────────────────────────────
# Times iEPPA vs L-BFGS-B at one (log_complexity, log_tol) point.
# Returns log(median_t_iEPPA / median_t_LBFGSB):
#   positive → iEPPA is slower → L-BFGS-B wins.
# Threshold: log(1.2) ≈ 0.182.
#
# n derivation: cats_per_margin fixed by complexity tercile, n = round(complexity/(K*cats)).
# Both solvers run with max_weight=Inf (unconstrained, exponential link).
# seed_extra: K-specific seed offset for K-stability sweeps (default 0L = main sweep)
#   K=3 stability: pass seed_extra = 3L * 10000000L
#   K=18 stability: pass seed_extra = 18L * 10000000L
time_cell <- function(log_complexity, log_tol, K = 9L, seed_extra = 0L) {
  K <- as.integer(K)
  cats_per_margin <- if (log_complexity <= 5.5) 4L else if (log_complexity <= 6.5) 8L else 16L
  n <- max(50L, 2L * K * cats_per_margin,
         as.integer(round(10^log_complexity / (K * cats_per_margin))))

  # Round-trip sanity check
  actual_log_c <- log10(n * K * cats_per_margin)
  if (abs(actual_log_c - log_complexity) > 0.2)
    warning(sprintf("time_cell: complexity round-trip %.2f log-units (lc=%.2f K=%d cats=%d n=%d)",
                    abs(actual_log_c - log_complexity), log_complexity, K, cats_per_margin, n))

  # Generate data; retry with progressively larger n if any margin has empty cells.
  # Empty cells make harvest() infeasible regardless of solver.
  seed_base <- bench_seed(log_complexity, log_tol) + seed_extra
  bd <- NULL
  n_try <- n
  for (.attempt in seq_len(5L)) {
    set.seed(seed_base + .attempt - 1L)
    cand <- make_bench_data(n_try, K, cats_per_margin)
    all_filled <- all(sapply(cand$df, function(col) min(table(col)) > 0L))
    if (all_filled) { bd <- cand; break }
    n_try <- as.integer(ceiling(n_try * 1.5))
  }
  if (is.null(bd)) stop("time_cell: could not generate feasible data after 5 attempts")
  conv <- list(absolute = 10^log_tol)

  time_algo <- function(method) {
    # 2 warmup runs (discarded)
    for (i in seq_len(2L))
      suppressWarnings(invisible(leafblower::harvest(
        bd$df, bd$targets, method = method,
        max_weight = Inf, min_weight = 0, convergence = conv, max_iterations = 500L)))
    # 5 timed runs
    median(replicate(5L, {
      t0 <- proc.time()[["elapsed"]]
      suppressWarnings(invisible(leafblower::harvest(
        bd$df, bd$targets, method = method,
        max_weight = Inf, min_weight = 0, convergence = conv, max_iterations = 500L)))
      proc.time()[["elapsed"]] - t0
    }))
  }

  t_ieppa  <- time_algo("ieppa")
  t_lbfgsb <- time_algo("lbfgsb")
  # Floor at 0.1ms to guard against proc.time() resolution yielding exact zeros.
  t_min <- 1e-4
  log(max(t_ieppa, t_min) / max(t_lbfgsb, t_min))
}

# ── fit_gp ────────────────────────────────────────────────────────────────────
# Fits a Matérn-5/2 GP to timing observations.
# design_mat: n_obs × 2 matrix (col 1 = log_complexity, col 2 = log_tol).
# y: numeric vector of log(t_iEPPA / t_LBFGSB).
fit_gp <- function(design_mat, y) {
  stopifnot(is.matrix(design_mat), nrow(design_mat) == length(y), ncol(design_mat) == 2L)
  DiceKriging::km(
    formula      = ~1,
    design       = as.data.frame(design_mat),
    response     = y,
    covtype      = "matern5_2",
    nugget.estim = TRUE,
    nugget       = 1e-4,          # lower bound: prevents degenerate fit on small n
    control      = list(trace = FALSE)
  )
}

# ── straddle_next ─────────────────────────────────────────────────────────────
# Straddle acquisition (Bryan et al. 2005): picks the candidate maximising
#   a(x) = -|μ(x) − threshold| + κ·σ(x)
# Pulls samples toward the contour (low |μ − threshold|) and uncertain regions (high σ).
straddle_next <- function(gp_model, candidates, threshold, kappa = 2) {
  pred <- DiceKriging::predict(gp_model,
                               newdata    = as.data.frame(candidates),
                               type       = "UK",
                               checkNames = FALSE)
  a   <- -abs(pred$mean - threshold) + kappa * pred$sd
  candidates[which.max(a), , drop = FALSE]
}

# ── classified_fraction ───────────────────────────────────────────────────────
# Fraction of candidates classified with ≥conf confidence as above or below threshold.
# Termination fires when this reaches 0.90.
classified_fraction <- function(gp_model, candidates, threshold, conf = 0.95) {
  pred    <- DiceKriging::predict(gp_model,
                                  newdata    = as.data.frame(candidates),
                                  type       = "UK",
                                  checkNames = FALSE)
  p_above <- pnorm(threshold, mean = pred$mean, sd = pred$sd, lower.tail = FALSE)
  mean(p_above > conf | p_above < (1 - conf))
}

# ── save_checkpoint ───────────────────────────────────────────────────────────
# Atomic checkpoint: write to .tmp, then rename to final path.
# state fields: design (n×2 matrix), y (numeric), gp (km or NULL),
#               iter (integer), classified (numeric), bounds (list).
save_checkpoint <- function(state, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(state, tmp)
  if (!file.rename(tmp, path))
    stop(sprintf("save_checkpoint: rename failed: %s -> %s", tmp, path))
  invisible(path)
}

# ── load_checkpoint ───────────────────────────────────────────────────────────
# Returns NULL if no checkpoint found; otherwise returns the saved state.
load_checkpoint <- function(path) {
  if (!file.exists(path)) return(NULL)
  readRDS(path)
}

# ── run_k_stability ───────────────────────────────────────────────────────────
# Evaluates a 4×4 fixed grid at K ∈ K_vals to check whether the 1.2× contour
# shifts with margin count. Overlays all contours on a single PDF.
#
# Decision rule: if the K=3 or K=18 contour shifts >0.5 log-units from the
# K=9 GP contour at any tol level, prints a scope warning.
run_k_stability <- function(state, K_vals = c(3L, 18L), threshold, out_dir = "benchmarks",
                             grid_size = 4L) {
  # grid_size: number of points per axis (default 4 → 4×4 production grid).
  #   Use grid_size=1 in tests to call a single time_cell() per K value.
  x1_pts <- seq(4.5, 7.2, length.out = grid_size)
  x2_pts <- seq(-5.5, -3.5, length.out = grid_size)
  grid   <- expand.grid(log_complexity = x1_pts, log_tol = x2_pts)

  results_list <- lapply(K_vals, function(K) {
    cat(sprintf("  K-stability: evaluating K=%d (16 points)...\n", K))
    y_K <- numeric(nrow(grid))
    for (i in seq_len(nrow(grid))) {
      # seed_extra = K * 10000000L ensures K-stability data is independent of main K=9 sweep.
      # The offset fits in 32-bit integer for K <= 214.
      # NOTE: seed_extra is passed INTO time_cell() — do NOT call set.seed() here.
      #       time_cell() internally calls set.seed(bench_seed(...) + seed_extra).
      y_K[i] <- tryCatch(
        time_cell(grid$log_complexity[i], grid$log_tol[i], K = K,
                  seed_extra = K * 10000000L),
        error = function(e) NA_real_
      )
    }
    data.frame(grid, y = y_K, K = K)
  })

  # Add K=9 GP posterior mean at the same grid points for comparison
  cands_grid <- as.matrix(grid)
  pred_k9    <- DiceKriging::predict(state$gp,
                                     newdata    = as.data.frame(cands_grid),
                                     type       = "UK",
                                     checkNames = FALSE)
  k9_df <- data.frame(grid, y = pred_k9$mean, K = 9L)

  all_df <- rbind(k9_df, do.call(rbind, results_list))
  all_df$K_label <- paste0("K=", all_df$K)

  # Contour shift check: compare each K_val contour against K=9 at each x2 level
  k9_vals <- k9_df$y
  for (K in K_vals) {
    kv <- all_df[all_df$K == K, "y"]
    if (any(is.finite(kv)) && any(is.finite(k9_vals))) {
      max_shift <- max(abs(kv[is.finite(kv) & is.finite(k9_vals)] -
                           k9_vals[is.finite(kv) & is.finite(k9_vals)]))
      if (max_shift > 0.5) {
        warning(sprintf(paste0(
          "K-stability: K=%d contour shifts %.2f log-units from K=9.\n",
          "Constants are valid only for K≈9 (Stepstone regime).\n",
          "Add comment to kComplexityThreshold/kTolThreshold in src/c_api.cpp:\n",
          "  // Threshold calibrated for K≈9, uniform-category problems."),
          K, max_shift))
      }
    }
  }

  # K-stability overlay plot
  p <- ggplot2::ggplot(all_df[is.finite(all_df$y), ],
                       ggplot2::aes(log_complexity, log_tol, z = y, colour = K_label)) +
    ggplot2::geom_contour(breaks = threshold, linewidth = 1) +
    ggplot2::labs(title = "K-stability: 1.2× contour at K=3, 9, 18",
                  x = "log10(complexity)", y = "log10(tol_abs)", colour = "K") +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(out_dir, "algo_selection_k_stability.pdf"), p, width = 8, height = 6)
  cat(sprintf("  K-stability plot saved to %s/algo_selection_k_stability.pdf\n", out_dir))
  invisible(all_df)
}

# ── run_benchmark ─────────────────────────────────────────────────────────────
# Main entry point. Runs the Bayesian LSE loop.
#
# budget:          max adaptive acquisitions (default 25)
# checkpoint_path: where to save/load state (default benchmarks/algo_selection_results.rds)
# out_dir:         where to write PDFs (default benchmarks/)
# seed:            seed for the 8-pt LHC initial design (default 42L)
run_benchmark <- function(budget          = 25L,
                          checkpoint_path = "benchmarks/algo_selection_results.rds",
                          out_dir         = "benchmarks",
                          seed            = 42L,
                          lhc_x1_max      = NULL) {
  # lhc_x1_max: cap LHC x1 (log_complexity) coordinates to at most this value.
  #   NULL (default) → full [4, 7.7] range used in production.
  #   5.5 → smoke-test mode: all LHC points stay in the fast sub-region
  #          (n ≤ ~4K at K=9, cats=8), keeping the smoke test to ~2 min.

  threshold  <- BENCH_THRESHOLD
  x1_range   <- BENCH_X1_RANGE
  x2_range   <- BENCH_X2_RANGE

  # 50×50 candidate grid — used for Straddle and classification
  candidates <- as.matrix(expand.grid(
    V1 = seq(x1_range[1], x1_range[2], length.out = 50L),
    V2 = seq(x2_range[1], x2_range[2], length.out = 50L)
  ))

  # ── Restart or initialise ──────────────────────────────────────────────────
  state <- load_checkpoint(checkpoint_path)
  if (!is.null(state)) {
    cat(sprintf("Restarting from checkpoint: %d evaluations, iter=%d, classified=%.2f\n",
                nrow(state$design), state$iter, state$classified))
    # Refit GP from saved design + y (do not re-run evaluations)
    state$gp <- fit_gp(state$design, state$y)
  } else {
    cat("Initialising: 8-point Latin hypercube design...\n")
    set.seed(seed)
    lhc_unit <- lhs::randomLHS(n = 8L, k = 2L)
    lhc_pts  <- cbind(
      x1_range[1] + lhc_unit[, 1L] * diff(x1_range),
      x2_range[1] + lhc_unit[, 2L] * diff(x2_range)
    )
    # Smoke-test cap: clamp LHC log_complexity coordinates to lhc_x1_max.
    # This keeps all initial evaluations in the fast sub-region without
    # changing the Straddle/classification grid (always full [4,7.7]).
    if (!is.null(lhc_x1_max))
      lhc_pts[, 1L] <- pmin(lhc_pts[, 1L], lhc_x1_max)
    design <- matrix(nrow = 0L, ncol = 2L)
    y      <- numeric(0L)
    for (i in seq_len(nrow(lhc_pts))) {
      cat(sprintf("  LHC %d/8: lc=%.2f, lt=%.2f ... ", i, lhc_pts[i, 1L], lhc_pts[i, 2L]))
      yi <- tryCatch(time_cell(lhc_pts[i, 1L], lhc_pts[i, 2L]),
                     error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NA_real_ })
      cat(sprintf("y=%s\n", if (is.finite(yi)) sprintf("%.3f", yi) else "FAILED"))
      design <- rbind(design, lhc_pts[i, , drop = FALSE])
      y      <- c(y, yi)
    }
    # Drop any NA evaluations
    ok     <- is.finite(y)
    design <- design[ok, , drop = FALSE]
    y      <- y[ok]

    if (sum(ok) < 2L)
      stop("All LHC evaluations failed — check leafblower installation and system capacity. ",
           "Cannot fit GP on fewer than 2 observations.")

    state <- list(design = design, y = y, gp = fit_gp(design, y),
                  iter = 0L, classified = 0.0, bounds = list(x1 = x1_range, x2 = x2_range))
    save_checkpoint(state, checkpoint_path)
  }

  # ── Adaptive acquisitions ──────────────────────────────────────────────────
  remaining <- max(0L, as.integer(budget) - state$iter)
  for (i in seq_len(remaining)) {
    state$classified <- classified_fraction(state$gp, candidates, threshold)
    cat(sprintf("Iter %d/%d: classified=%.2f\n", state$iter + 1L, budget, state$classified))
    if (state$classified >= 0.90) {
      cat("Termination: 90% classified.\n")
      state$converged <- TRUE
      break
    }

    next_pt <- straddle_next(state$gp, candidates, threshold)
    cat(sprintf("  Next: lc=%.3f, lt=%.3f ... ", next_pt[1L, 1L], next_pt[1L, 2L]))
    yi <- tryCatch(time_cell(next_pt[1L, 1L], next_pt[1L, 2L]),
                   error = function(e) { cat("ERROR:", conditionMessage(e), "\n"); NA_real_ })
    cat(sprintf("y=%s\n", if (is.finite(yi)) sprintf("%.3f", yi) else "FAILED"))

    if (is.finite(yi)) {
      state$design <- rbind(state$design, next_pt)
      state$y      <- c(state$y, yi)
      state$gp     <- fit_gp(state$design, state$y)
    }
    state$iter <- state$iter + 1L

    if (state$iter %% 5L == 0L)
      save_checkpoint(state, checkpoint_path)
  }

  # Final checkpoint
  save_checkpoint(state, checkpoint_path)

  # ── Poor-fit warning ───────────────────────────────────────────────────────
  final_classified <- classified_fraction(state$gp, candidates, threshold)
  if (final_classified < 0.90) {
    warning(
      sprintf(paste0(
        "GP classification %.0f%% < 90%% at termination.\n",
        "Inspect benchmarks/algo_selection_uncertainty.pdf before committing constants.\n",
        "Consider filing a follow-up issue:\n",
        "  bd create --title='algo-selection: 3D sweep needed (K/n confound)' ",
        "--description='GP classified only %.0f%% of space after %d acquisitions' ",
        "--type=task --priority=3"),
        100 * final_classified, 100 * final_classified, state$iter))
  }

  # ── Plots ─────────────────────────────────────────────────────────────────
  cat("Generating plots...\n")
  make_plots(state, candidates, threshold, out_dir)

  cat(sprintf("\nDone. %d total evaluations, %.0f%% classified.\n",
              nrow(state$design), 100 * final_classified))
  cat(sprintf("Inspect: %s/algo_selection_contour.pdf\n", out_dir))
  invisible(state)
}

# ── Entry point ───────────────────────────────────────────────────────────────
# Runs only when executed via `Rscript benchmarks/algo_selection_benchmark.R`.
# Skipped when sourced in tests (set `.BENCH_SOURCED <- TRUE` before source()).
# K-stability is called here (not inside run_benchmark) so smoke tests can call
# run_benchmark() alone without triggering 32 extra time_cell() evaluations.
if (!.BENCH_SOURCED) {
  state <- run_benchmark()
  run_k_stability(state, K_vals = c(3L, 18L), threshold = BENCH_THRESHOLD,
                  out_dir = "benchmarks")
}
