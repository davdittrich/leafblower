#!/usr/bin/env Rscript
# benchmarks/parity_compare.R
# Cross-language parity report.
# Reads R + Python metric tables, joins by method, computes per-method
# weight-vector diff statistics, and prints a single table sorted by R wall_time.

suppressPackageStartupMessages({ library(arrow) })

R_CSV   <- "benchmarks/results/parity/parity_bench_R.csv"
PY_CSV  <- "benchmarks/results/parity/parity_bench_Py.csv"
WDIR_R  <- "benchmarks/results/parity/weights_R"
WDIR_PY <- "benchmarks/results/parity/weights_Py"
OUT_CSV <- "benchmarks/results/parity/parity_bench_combined.csv"

slug <- function(s) gsub("[^A-Za-z0-9]+", "_", s)

if (!file.exists(R_CSV))  stop("Missing ", R_CSV,  " — run parity_bench.R first")
if (!file.exists(PY_CSV)) stop("Missing ", PY_CSV, " — run parity_bench.py first")

R  <- read.csv(R_CSV,  stringsAsFactors = FALSE)
Py <- read.csv(PY_CSV, stringsAsFactors = FALSE)

methods <- intersect(R$method, Py$method)
miss_R  <- setdiff(Py$method, R$method)
miss_Py <- setdiff(R$method,  Py$method)
if (length(miss_R) > 0)  cat("WARN: missing in R:  ", paste(miss_R,  collapse=", "), "\n")
if (length(miss_Py) > 0) cat("WARN: missing in Py: ", paste(miss_Py, collapse=", "), "\n")

diff_stats <- function(label) {
  fp_R  <- file.path(WDIR_R,  paste0(slug(label), ".feather"))
  fp_Py <- file.path(WDIR_PY, paste0(slug(label), ".feather"))
  if (!file.exists(fp_R) || !file.exists(fp_Py)) {
    return(c(max_abs_diff = NA_real_, mean_abs_diff = NA_real_,
             rel_max_diff = NA_real_, identical_1e4 = NA, n_R = NA, n_Py = NA))
  }
  wR  <- as.numeric(read_feather(fp_R)$w)
  wPy <- as.numeric(read_feather(fp_Py)$w)
  if (length(wR) != length(wPy)) {
    return(c(max_abs_diff = NA_real_, mean_abs_diff = NA_real_,
             rel_max_diff = NA_real_, identical_1e4 = FALSE,
             n_R = length(wR), n_Py = length(wPy)))
  }
  d <- abs(wR - wPy)
  denom <- pmax(abs(wR), 1e-12)
  c(max_abs_diff  = max(d),
    mean_abs_diff = mean(d),
    rel_max_diff  = max(d / denom),
    identical_1e4 = isTRUE(max(d) < 1e-4),
    n_R = length(wR), n_Py = length(wPy))
}

cat("\nComparing weight vectors per method...\n")
diffs <- t(sapply(methods, diff_stats))
diffs <- as.data.frame(diffs, stringsAsFactors = FALSE)
diffs$method <- methods

# Join everything on method
colnames(R)[-1]  <- paste0(colnames(R)[-1],  "_R")
colnames(Py)[-1] <- paste0(colnames(Py)[-1], "_Py")
J <- merge(R, Py, by = "method", all = TRUE)
J <- merge(J, diffs, by = "method", all = TRUE)

# Sort by R wall time (median_ms_R), then Py if R missing
J$sort_key <- ifelse(!is.na(J$median_ms_R), J$median_ms_R, J$median_ms_Py)
J <- J[order(J$sort_key), , drop = FALSE]
J$sort_key <- NULL

# Identical-flag readability
J$identical_1e4 <- ifelse(is.na(J$identical_1e4), "NA",
                   ifelse(as.logical(J$identical_1e4), "YES", "NO"))

write.csv(J, OUT_CSV, row.names = FALSE)
cat(sprintf("Saved combined table: %s\n\n", OUT_CSV))

# --- Console report ----------------------------------------------------------
cat("=== Parity Report — sorted by R median_ms ===\n\n")

# Wall-time + status table
wt <- data.frame(
  method      = J$method,
  R_ms        = J$median_ms_R,
  Py_ms       = J$median_ms_Py,
  ratio_PyR   = round(J$median_ms_Py / pmax(J$median_ms_R, 1), 2),
  R_iter      = J$iterations_R,
  Py_iter     = J$iterations_Py,
  R_status    = J$status_R,
  Py_status   = J$status_Py,
  stringsAsFactors = FALSE
)
cat("--- Wall time, iterations, status ---\n")
print(wt, row.names = FALSE)

# Metrics table
mt <- data.frame(
  method      = J$method,
  R_max_err   = signif(J$max_err_R, 4),
  Py_max_err  = signif(J$max_err_Py, 4),
  R_mkl       = signif(J$marginal_kl_R, 4),
  Py_mkl      = signif(J$marginal_kl_Py, 4),
  R_kl        = signif(J$kl_R, 4),
  Py_kl       = signif(J$kl_Py, 4),
  R_chi2      = signif(J$chi2_R, 4),
  Py_chi2     = signif(J$chi2_Py, 4),
  stringsAsFactors = FALSE
)
cat("\n--- All metrics ---\n")
print(mt, row.names = FALSE)

# Weight parity
wp <- data.frame(
  method        = J$method,
  max_abs_diff  = signif(as.numeric(J$max_abs_diff), 4),
  mean_abs_diff = signif(as.numeric(J$mean_abs_diff), 4),
  rel_max_diff  = signif(as.numeric(J$rel_max_diff), 4),
  identical_1e4 = J$identical_1e4,
  stringsAsFactors = FALSE
)
cat("\n--- Weight-vector R↔Py diff (tol=1e-4) ---\n")
print(wp, row.names = FALSE)

n_yes <- sum(wp$identical_1e4 == "YES", na.rm = TRUE)
n_no  <- sum(wp$identical_1e4 == "NO",  na.rm = TRUE)
n_na  <- sum(wp$identical_1e4 == "NA",  na.rm = TRUE)
cat(sprintf("\nVerdict: %d/%d methods identical at tol=1e-4 (%d differ, %d unavailable)\n",
            n_yes, length(methods), n_no, n_na))
