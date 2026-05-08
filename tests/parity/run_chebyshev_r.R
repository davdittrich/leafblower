#!/usr/bin/env Rscript
# run_chebyshev_r.R — T5 parity helper: chebyshev at tol=1e-4, max_weight=5
# Args: <data_csv> <targets_json> <out_csv>
# out_csv columns: weight, iterations, status

# Enforce single-thread BLAS before any library load so that thread count is
# locked before OpenBLAS / MKL / OpenMP initialise their thread pools.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1")

tryCatch(
  suppressPackageStartupMessages({
    library(leafblower)
    library(jsonlite)
  }),
  error = function(e) {
    cat("SKIP:", conditionMessage(e), "\n", file = stderr())
    quit(status = 2, save = "no")
  }
)

args <- commandArgs(trailingOnly = TRUE)
df   <- read.csv(args[1], stringsAsFactors = FALSE)
tgt  <- fromJSON(args[2])
# Normalize each variable's proportions (JSON may not be normalized)
tgt  <- lapply(tgt, function(x) { v <- unlist(x); v / sum(v) })

res <- harvest(
  df, tgt,
  method         = "chebyshev",
  min_weight     = 0,
  max_weight     = 5,
  max_iterations = 3000,
  convergence    = list(metric = "max_err", rule = "improvement", tol = 1e-4),
  verbose        = 0L
)

w      <- res[["weights"]]
r      <- attr(res, "result")
iters  <- r[["iterations"]]
status <- r[["status"]]

write.csv(
  data.frame(weight = w, iterations = iters, status = status),
  args[3], row.names = FALSE
)
