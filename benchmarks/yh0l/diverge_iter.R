#!/usr/bin/env Rscript
# benchmarks/yh0l/diverge_iter.R
# T2: load trajectory.csv; locate first iter where |convergence_metric_R - convergence_metric_Py| > 1e-12.
# Also parses terminal status + iter from log files.

suppressPackageStartupMessages({
  if (!requireNamespace("readr", quietly = TRUE)) {
    read_csv <- function(path, ...) utils::read.csv(path, stringsAsFactors = FALSE)
  } else {
    library(readr)
  }
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args, value = TRUE)
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_flag[1])))
} else {
  SCRIPT_DIR <- getwd()
}
ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)

CSV    <- file.path(ROOT, "benchmarks", "yh0l", "trajectory.csv")
LOG_R  <- file.path(ROOT, "benchmarks", "yh0l", "traj_R.log")
LOG_PY <- file.path(ROOT, "benchmarks", "yh0l", "traj_Py.log")

# ── Parse termination line from log ─────────────────────────────────────────
parse_terminate <- function(path) {
  if (!file.exists(path)) return(list(status = "LOG_MISSING", iterations = NA_integer_, max_error = NA_real_))
  lines <- readLines(path, warn = FALSE)

  # Primary: explicit TERMINATE line written by run scripts
  pat <- "## TERMINATE status=(\\S+) \\((\\S+)\\) iterations=(\\S+) max_error=(\\S+)"
  for (ln in rev(lines)) {
    m <- regexec(pat, ln, perl = TRUE)
    parts <- regmatches(ln, m)[[1]]
    if (length(parts) > 1) {
      return(list(
        status     = parts[3],
        iterations = suppressWarnings(as.integer(parts[4])),
        max_error  = suppressWarnings(as.numeric(parts[5]))
      ))
    }
  }

  # Parse TERMINATE line: "## TERMINATE status=<s> (<label>) iterations=<n> max_error=<v>"
  conv_pat  <- "## TERMINATE status=\\S+ \\(\\S+\\) iterations=(\\d+|NA) max_error=([0-9eE.+\\-]+|NA)"
  for (ln in rev(lines)) {
    m <- regexec(conv_pat, ln, perl = TRUE)
    parts <- regmatches(ln, m)[[1]]
    if (length(parts) > 1) {
      verb  <- parts[2]
      iters <- suppressWarnings(as.integer(parts[3]))
      err   <- suppressWarnings(as.numeric(parts[4]))
      st    <- if (grepl("converge", verb)) "converged" else "no_conv"
      return(list(status = st, iterations = iters, max_error = err))
    }
  }

  list(status = "NOT_FOUND", iterations = NA_integer_, max_error = NA_real_)
}

term_R  <- parse_terminate(LOG_R)
term_Py <- parse_terminate(LOG_PY)
cat(sprintf("[diverge_iter] R  terminate: status=%s iter=%s max_err=%s\n",
            term_R$status, term_R$iterations, term_R$max_error))
cat(sprintf("[diverge_iter] Py terminate: status=%s iter=%s max_err=%s\n",
            term_Py$status, term_Py$iterations, term_Py$max_error))

# ── Load CSV ─────────────────────────────────────────────────────────────────
if (!file.exists(CSV)) stop(sprintf("trajectory.csv not found: %s", CSV))
traj <- read_csv(CSV, show_col_types = FALSE)
cat(sprintf("[diverge_iter] CSV rows: %d  (R:%d  Py:%d)\n",
            nrow(traj), sum(traj$lang == "R"), sum(traj$lang == "Py")))

# ── Align on common iters ────────────────────────────────────────────────────
df_r  <- traj[traj$lang == "R",  c("iter", "convergence_metric")]
df_py <- traj[traj$lang == "Py", c("iter", "convergence_metric")]
joined <- merge(df_r, df_py, by = "iter", suffixes = c("_R", "_Py"))
joined <- joined[order(joined$iter), ]

if (nrow(joined) == 0) {
  cat("[diverge_iter] No overlapping iters — cannot compute divergence.\n")
  cat("diverge_iter: NULL\n")
  quit(save = "no")
}

joined$abs_diff <- abs(joined$convergence_metric_R - joined$convergence_metric_Py)
joined$ratio    <- joined$convergence_metric_R / joined$convergence_metric_Py

THRESHOLD <- 1e-12
diverged  <- joined[joined$abs_diff > THRESHOLD, ]

if (nrow(diverged) == 0) {
  cat(sprintf("[diverge_iter] No divergence > %.0e across %d shared iters.\n",
              THRESHOLD, nrow(joined)))
  cat(sprintf("diverge_iter: NULL  (all |diff| <= %.0e)\n", THRESHOLD))
  cat(sprintf("max_abs_diff: %.6e at iter %d\n",
              max(joined$abs_diff), joined$iter[which.max(joined$abs_diff)]))
} else {
  first <- diverged[1, ]
  cat(sprintf("\n[diverge_iter] First divergence > 1e-12 at iter %d:\n", first$iter))
  cat(sprintf("  convergence_metric_R  = %.15e\n", first$convergence_metric_R))
  cat(sprintf("  convergence_metric_Py = %.15e\n", first$convergence_metric_Py))
  cat(sprintf("  abs_diff = %.6e\n",  first$abs_diff))
  cat(sprintf("  ratio    = %.10f\n", first$ratio))
  cat(sprintf("\ndivergence_iter: %d\n", first$iter))
  cat(sprintf("metric_R_at_diverge:  %.15e\n", first$convergence_metric_R))
  cat(sprintf("metric_Py_at_diverge: %.15e\n", first$convergence_metric_Py))
}

cat(sprintf("\nstatus_R_at_terminate:  %s\n", term_R$status))
cat(sprintf("status_Py_at_terminate: %s\n", term_Py$status))
cat(sprintf("iter_R_terminate:  %s\n", term_R$iterations))
cat(sprintf("iter_Py_terminate: %s\n", term_Py$iterations))

# First 10 diverged rows
if (nrow(diverged) > 0) {
  cat("\n[diverge_iter] First 10 diverged iters:\n")
  print(head(diverged[, c("iter","convergence_metric_R","convergence_metric_Py","abs_diff","ratio")], 10))
}
