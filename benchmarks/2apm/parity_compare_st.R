#!/usr/bin/env Rscript
# benchmarks/2apm/parity_compare_st.R
# T2 (leafblower-2apm.2): compare single-thread R vs Py chebyshev weights.
# Reads 2apm/ CSVs and weight feathers, computes d_st, writes parity_singlethread.csv.

suppressPackageStartupMessages({ library(arrow) })

R_CSV   <- "benchmarks/2apm/parity_bench_st_R.csv"
PY_CSV  <- "benchmarks/2apm/parity_bench_st_Py.csv"
WDIR_R  <- "benchmarks/2apm/weights_R_st"
WDIR_PY <- "benchmarks/2apm/weights_Py_st"
OUT_CSV <- "benchmarks/2apm/parity_singlethread.csv"

slug <- function(s) gsub("[^A-Za-z0-9]+", "_", s)

if (!file.exists(R_CSV))  stop("Missing ", R_CSV,  " — run parity_bench_st_R.R first")
if (!file.exists(PY_CSV)) stop("Missing ", PY_CSV, " — run parity_bench_st_py.py first")

R  <- read.csv(R_CSV,  stringsAsFactors = FALSE)
Py <- read.csv(PY_CSV, stringsAsFactors = FALSE)

methods <- intersect(R$method, Py$method)

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

cat("\nComputing single-thread weight diffs...\n")
diffs <- t(sapply(methods, diff_stats))
diffs <- as.data.frame(diffs, stringsAsFactors = FALSE)
diffs$method <- methods

R_sub  <- R[R$method  %in% methods, ]
Py_sub <- Py[Py$method %in% methods, ]
names(R_sub)  <- paste0(names(R_sub),  "_R")
names(Py_sub) <- paste0(names(Py_sub), "_Py")
R_sub$method  <- methods
Py_sub$method <- methods

J <- merge(R_sub,  diffs,  by = "method")
J <- merge(J,      Py_sub, by = "method")
J$thread_mode <- "single"

# Bring in multi-thread d_mt from existing combined CSV for chebyshev
MT_CSV <- "benchmarks/results/parity/parity_bench_combined.csv"
if (file.exists(MT_CSV)) {
  mt <- read.csv(MT_CSV, stringsAsFactors = FALSE)
  cheb_mt <- mt[mt$method == "chebyshev", "max_abs_diff"]
  if (length(cheb_mt) == 1) {
    J$d_mt_chebyshev <- cheb_mt
  }
}

J$identical_1e4 <- ifelse(is.na(J$identical_1e4), "NA",
                   ifelse(as.logical(J$identical_1e4), "YES", "NO"))

write.csv(J, OUT_CSV, row.names = FALSE)
cat(sprintf("Saved: %s\n", OUT_CSV))

# Summary
cat("\n=== T2 Cause Attribution: Chebyshev ===\n")
cheb <- J[J$method == "chebyshev", ]
d_st <- as.numeric(cheb$max_abs_diff)
d_mt <- if ("d_mt_chebyshev" %in% names(cheb)) as.numeric(cheb$d_mt_chebyshev) else NA

cat(sprintf("  d_st (single-thread max_abs_diff): %.4e\n", d_st))
cat(sprintf("  d_mt (multi-thread  max_abs_diff): %.4e\n", d_mt))

if (!is.na(d_st)) {
  if (d_st < 1e-9 && !is.na(d_mt) && d_mt > 1e-3) {
    cat("  => cause (d): BLAS thread nondeterminism\n")
  } else if (d_st > 1e-3) {
    cat("  => cause (c): shared-code numerical drift\n")
  } else {
    cat(sprintf("  => cause (e): AMBIGUOUS (1e-9 <= d_st=%.4e <= 1e-3)\n", d_st))
    cat("  => Run NREPS=5 sub-resolution (see parity_compare_st_nreps.R)\n")
  }
}
