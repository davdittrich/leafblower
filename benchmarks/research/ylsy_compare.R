# WU-7: Unified comparison harness for ylsy CP+IPM spike investigation.
#
# Reads CP and IPM summary CSVs (WU-3, WU-6) and regenerates baselines INLINE
# via harvest() (ieppa+sraa, newton_kl, lbfgsb) on each fixture, producing
# one unified comparison CSV. Also emits log-log trajectory plots and computes
# the R3 falsification result on kk1204.
#
# Run from package root:  Rscript benchmarks/research/ylsy_compare.R

suppressPackageStartupMessages({
  library(leafblower)
  library(arrow)
  library(jsonlite)
  library(RhpcBLASctl)
})

Sys.setenv(OMP_NUM_THREADS = "1")
blas_set_num_threads(1L)

source("benchmarks/research/utils.R")

stopifnot(file.exists("benchmarks/stepstone_bench_data.parquet"))
stopifnot(file.exists("benchmarks/stepstone_bench_targets.json"))

git_sha <- system("git rev-parse --short HEAD", intern = TRUE)
omp_threads <- Sys.getenv("OMP_NUM_THREADS")
blas_threads <- blas_get_num_procs()

results_dir <- "benchmarks/research/results"
plots_dir <- file.path(results_dir, "plots")
dir.create(plots_dir, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# ---------------------------------------------------------------------------
# Read CP + IPM summaries (WU-3, WU-6)
# ---------------------------------------------------------------------------
cp_summary <- read.csv(file.path(results_dir, "cp_summary.csv"), stringsAsFactors = FALSE)
ipm_summary <- read.csv(file.path(results_dir, "ipm_summary.csv"), stringsAsFactors = FALSE)

# ---------------------------------------------------------------------------
# Hermetic baselines via harvest(): ieppa (accelerate), newton_kl, lbfgsb
# ---------------------------------------------------------------------------
schema_cols <- colnames(cp_summary)

baseline_row <- function(solver, fixture, rep, res, wall_s, fx) {
  status_code <- if (is.null(res$status)) 1L else as.integer(res$status)
  max_err <- if (is.null(res$max_error)) NA_real_ else as.numeric(res$max_error)
  n_iter <- if (is.null(res$iterations)) NA_integer_ else as.integer(res$iterations)
  status_msg <- switch(as.character(status_code),
                       "0" = "converged",
                       "1" = "did not converge",
                       "2" = "infeasible",
                       "3" = "bad argument",
                       "4" = "budget exhausted",
                       "5" = "loss plateau",
                       "unknown")
  parity <- if (fixture == "stepstone_K9" && is.finite(max_err)) max_err / 1.13e-4 else NA_real_

  data.frame(
    solver = solver, fixture = fixture, rep = rep,
    status_code = status_code, status_msg = status_msg,
    converged = (status_code == 0L),
    max_err = max_err, max_err_ergodic = NA_real_,
    wall_s = wall_s, n_iter = n_iter,
    final_primal_resid = NA_real_, final_residual = NA_real_,
    n_projected_dims_max = NA_real_,
    rate_exponent_last = NA_real_, rate_R2_last = NA_real_, n_fit_points_last = NA_integer_,
    rate_exponent_ergodic = NA_real_, rate_R2_ergodic = NA_real_, n_fit_points_ergodic = NA_integer_,
    stepstone_parity_ratio = parity,
    seed = 1L, omp_threads = omp_threads, blas_threads = blas_threads, git_sha = git_sha,
    stringsAsFactors = FALSE
  )
}

run_baseline <- function(method, fixture, fx, rep) {
  t0 <- Sys.time()
  r <- tryCatch(
    suppressWarnings(harvest(
      fx$df, fx$tgt,
      method = method,
      max_weight = fx$max_weight,
      attach_weights = FALSE,
      max_iterations = if (method == "ieppa") 50L else 200L,
      verbose = 0L,
      accelerate = (method == "ieppa")
    )),
    error = function(e) NULL
  )
  wall_s <- as.numeric(Sys.time() - t0, units = "secs")
  if (is.null(r)) {
    res <- list(status = 3L, max_error = NA_real_, iterations = NA_integer_)
  } else {
    res <- attr(r, "result")
    if (is.null(res)) res <- list(status = 1L, max_error = NA_real_, iterations = NA_integer_)
  }
  solver_name <- if (method == "ieppa") "ieppa+sraa" else method
  baseline_row(solver_name, fixture, rep, res, wall_s, fx)
}

cat("WU-7: regenerating baselines via harvest()\n")
fixtures <- list(
  t1_small     = make_t1_small(seed = 1L),
  stepstone_K9 = make_stepstone_K9(),
  kk1204_K20   = make_kk1204_K20(seed = 1L)
)

baseline_rows <- list()
for (fx_name in names(fixtures)) {
  fx <- fixtures[[fx_name]]
  reps <- if (fx_name == "kk1204_K20") 3L else 1L
  for (method in c("ieppa", "newton_kl", "lbfgsb")) {
    for (rep in seq_len(reps)) {
      cat(sprintf("  baseline: %s on %s rep=%d ...\n", method, fx_name, rep))
      row <- run_baseline(method, fx_name, fx, rep)
      cat(sprintf("    -> max_err=%.3e wall=%.2fs status=%d\n",
                  row$max_err %||% NA, row$wall_s, row$status_code))
      baseline_rows[[length(baseline_rows) + 1L]] <- row
    }
  }
}

baseline_df <- do.call(rbind, baseline_rows)

# Reorder to match CP/IPM schema and rbind
unified <- rbind(cp_summary[, schema_cols], ipm_summary[, schema_cols], baseline_df[, schema_cols])
write.csv(unified, file.path(results_dir, "ylsy_comparison.csv"), row.names = FALSE)
cat(sprintf("WU-7: ylsy_comparison.csv written (%d rows)\n", nrow(unified)))

# ---------------------------------------------------------------------------
# Trajectory plots: log10(max_err) vs log10(iter) for CP and IPM per fixture
# ---------------------------------------------------------------------------
have_ggplot2 <- requireNamespace("ggplot2", quietly = TRUE)

plot_fixture <- function(fx_name) {
  cp_path <- file.path(results_dir, sprintf("cp_%s_trace_rep1.csv", fx_name))
  ipm_path <- file.path(results_dir, sprintf("ipm_%s_trace_rep1.csv", fx_name))
  cp_tr  <- if (file.exists(cp_path))  read.csv(cp_path,  stringsAsFactors = FALSE) else NULL
  ipm_tr <- if (file.exists(ipm_path)) read.csv(ipm_path, stringsAsFactors = FALSE) else NULL

  rows <- list()
  if (!is.null(cp_tr) && nrow(cp_tr) > 0) {
    keep <- cp_tr$iter > 0 & is.finite(cp_tr$max_err_last) & cp_tr$max_err_last > 0
    if (any(keep)) rows[[length(rows) + 1L]] <- data.frame(
      iter = cp_tr$iter[keep], max_err = cp_tr$max_err_last[keep],
      solver = "cp_last", stringsAsFactors = FALSE)
    keep_e <- cp_tr$iter > 0 & is.finite(cp_tr$max_err_ergodic) & cp_tr$max_err_ergodic > 0
    if (any(keep_e)) rows[[length(rows) + 1L]] <- data.frame(
      iter = cp_tr$iter[keep_e], max_err = cp_tr$max_err_ergodic[keep_e],
      solver = "cp_ergodic", stringsAsFactors = FALSE)
  }
  if (!is.null(ipm_tr) && nrow(ipm_tr) > 0) {
    iter <- if ("iter_outer" %in% names(ipm_tr)) ipm_tr$iter_outer + 1L else seq_len(nrow(ipm_tr))
    keep <- is.finite(ipm_tr$max_err) & ipm_tr$max_err > 0
    if (any(keep)) rows[[length(rows) + 1L]] <- data.frame(
      iter = iter[keep], max_err = ipm_tr$max_err[keep],
      solver = "ipm", stringsAsFactors = FALSE)
  }
  if (length(rows) == 0) return(invisible(NULL))
  long <- do.call(rbind, rows)

  out <- file.path(plots_dir, sprintf("%s_trajectory.png", fx_name))
  if (have_ggplot2) {
    p <- ggplot2::ggplot(long, ggplot2::aes(x = iter, y = max_err, colour = solver)) +
      ggplot2::geom_line() + ggplot2::geom_point(size = 0.6) +
      ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
      ggplot2::labs(x = "iteration (log)", y = "max_err (log)",
                    title = sprintf("%s: CP + IPM trajectory", fx_name),
                    colour = "solver") +
      ggplot2::theme_minimal()
    ggplot2::ggsave(out, p, width = 6, height = 4, dpi = 120)
  } else {
    png(out, width = 720, height = 480)
    plot(NA, xlim = range(long$iter), ylim = range(long$max_err),
         log = "xy", xlab = "iteration", ylab = "max_err",
         main = sprintf("%s: CP + IPM trajectory", fx_name))
    cols <- c(cp_last = "black", cp_ergodic = "blue", ipm = "red")
    for (s in unique(long$solver)) {
      sub <- long[long$solver == s, ]
      lines(sub$iter, sub$max_err, col = cols[[s]] %||% "gray")
      points(sub$iter, sub$max_err, col = cols[[s]] %||% "gray", pch = 20, cex = 0.4)
    }
    legend("topright", legend = unique(long$solver),
           col = cols[unique(long$solver)], lty = 1, bty = "n")
    dev.off()
  }
  cat(sprintf("WU-7: plot %s written\n", out))
  invisible(out)
}

for (fx_name in c("kk1204_K20", "stepstone_K9", "t1_small")) plot_fixture(fx_name)

# ---------------------------------------------------------------------------
# R3 falsification: 1% agreement across CP, IPM, ieppa+sraa, newton_kl on kk1204
# ---------------------------------------------------------------------------
kk_rows <- unified[unified$fixture == "kk1204_K20" & unified$rep == 1L, ]
r3_solvers <- c("cp", "ipm", "ieppa+sraa", "newton_kl")
r3_vals <- setNames(rep(NA_real_, length(r3_solvers)), r3_solvers)
for (s in r3_solvers) {
  v <- kk_rows$max_err[kk_rows$solver == s]
  if (length(v) >= 1L) r3_vals[s] <- v[1L]
}
finite_vals <- r3_vals[is.finite(r3_vals)]
if (length(finite_vals) >= 2L) {
  ratio <- max(finite_vals) / min(finite_vals)
  r3_status <- if (ratio <= 1.01) "CONFIRMED" else "RULED_OUT"
} else {
  ratio <- NA_real_
  r3_status <- "INSUFFICIENT_DATA"
}
cat(sprintf("WU-7: R3 falsification: %s (max/min = %.3f); values:\n", r3_status, ratio))
print(r3_vals)

writeLines(c(
  sprintf("status=%s", r3_status),
  sprintf("ratio=%.6f", ratio),
  sprintf("cp=%.6e", r3_vals["cp"]),
  sprintf("ipm=%.6e", r3_vals["ipm"]),
  sprintf("ieppa_sraa=%.6e", r3_vals["ieppa+sraa"]),
  sprintf("newton_kl=%.6e", r3_vals["newton_kl"])
), con = file.path(results_dir, "r3_falsification.txt"))

cat("WU-7 ylsy_compare.R complete\n")
