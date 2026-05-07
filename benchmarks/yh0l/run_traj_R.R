#!/usr/bin/env Rscript
# benchmarks/yh0l/run_traj_R.R
# T2: run ieppa_soft with verbose=2 on stepstone fulldata; save log to traj_R.log
# Run with: OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 Rscript benchmarks/yh0l/run_traj_R.R

suppressPackageStartupMessages({
  library(arrow)
  library(jsonlite)
  library(leafblower)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# Resolve script location robustly under Rscript
args <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args, value = TRUE)
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_flag[1])))
} else {
  SCRIPT_DIR <- getwd()
}
ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)

DATA_PATH    <- file.path(ROOT, "benchmarks", "stepstone_fulldata_bench_data.parquet")
TARGETS_PATH <- file.path(ROOT, "benchmarks", "stepstone_fulldata_bench_targets.json")
LOG_PATH     <- file.path(ROOT, "benchmarks", "yh0l", "traj_R.log")

cat(sprintf("[run_traj_R] ROOT=%s\n", ROOT))

# ── Load data ─────────────────────────────────────────────────────────────────
df_raw  <- as.data.frame(read_parquet(DATA_PATH))
tgt_raw <- jsonlite::fromJSON(TARGETS_PATH)
tgt     <- lapply(tgt_raw, function(x) { x <- unlist(x); x / sum(x) })
df      <- df_raw[, names(tgt), drop = FALSE]
for (v in names(tgt)) df[[v]] <- factor(df[[v]], levels = names(tgt[[v]]))
cat(sprintf("[run_traj_R] n=%s margins=%d\n", format(nrow(df), big.mark = ","), length(tgt)))
cat(sprintf("[run_traj_R] writing verbose=2 log to %s\n", LOG_PATH))

# ── Sink stdout+stderr to log; run solver ──────────────────────────────────
con <- file(LOG_PATH, open = "wt")
sink(con, type = "output")
sink(con, type = "message")

cat(sprintf("## run_traj_R.R start: %s\n", Sys.time()))
cat(sprintf("## OMP_NUM_THREADS=%s OPENBLAS_NUM_THREADS=%s\n",
            Sys.getenv("OMP_NUM_THREADS", "(unset)"),
            Sys.getenv("OPENBLAS_NUM_THREADS", "(unset)")))

status_map <- c("0" = "converged", "1" = "no_conv", "2" = "infeasible",
                "3" = "bad_arg",   "4" = "budget",   "5" = "stall")

res <- tryCatch(
  withCallingHandlers(
    leafblower::harvest(
      df, tgt,
      method         = "ieppa_soft",
      convergence    = list(tol = 1e-4),
      max_weight     = 5,
      max_iterations = 3000L,
      attach_weights = FALSE,
      verbose        = 2L
    ),
    warning = function(w) {
      cat(sprintf("## WARNING: %s\n", conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  ),
  error = function(e) {
    cat(sprintf("## ERROR: %s\n", conditionMessage(e)))
    NULL
  }
)

# attach_weights=FALSE → numeric vector with attr(res, "result") = calib_result
cr <- if (!is.null(res)) attr(res, "result") else NULL
if (!is.null(cr) && !is.null(cr$status)) {
  st <- as.character(cr$status)
  cat(sprintf("## TERMINATE status=%s (%s) iterations=%d max_error=%.6e\n",
              st, status_map[[st]] %||% "unknown", cr$iterations, cr$max_error))
} else {
  cat("## TERMINATE status=ERROR iterations=NA max_error=NA\n")
}
cat(sprintf("## run_traj_R.R end: %s\n", Sys.time()))

sink(type = "output")
sink(type = "message")
close(con)

cat(sprintf("[run_traj_R] done. Log: %s\n", LOG_PATH))
