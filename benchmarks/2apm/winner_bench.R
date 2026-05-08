#!/usr/bin/env Rscript
# benchmarks/2apm/winner_bench.R
# T3 (leafblower-2apm.3): multi-DGP chebyshev winner determination — R side.
# Run AFTER winner_bench.py (single-thread BLAS sequential to avoid oversubscription).
# Outputs: benchmarks/2apm/results_R_<dgp>.csv + weights feather per DGP/rep.

# ── Single-thread BLAS (must be before library load) ──────────────────────────
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")

suppressPackageStartupMessages({
  library(arrow)
  library(jsonlite)
  library(leafblower)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ── Locate repo root ──────────────────────────────────────────────────────────
args        <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args, value = TRUE)
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_flag[1])))
} else {
  SCRIPT_DIR <- getwd()
}
ROOT   <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)
APMDIR <- file.path(ROOT, "benchmarks", "2apm")
WDIR   <- file.path(APMDIR, "weights_R")
dir.create(WDIR, recursive = TRUE, showWarnings = FALSE)

NREPS  <- 3L
MAX_IT <- 3000L
MAX_WT <- 5.0
TOL    <- 1e-4

# ── DGP definitions ───────────────────────────────────────────────────────────
dgps <- list(
  fulldata = list(
    data_path    = file.path(ROOT, "benchmarks", "stepstone_fulldata_bench_data.parquet"),
    targets_path = file.path(ROOT, "benchmarks", "stepstone_fulldata_bench_targets.json")
  ),
  medium = list(
    data_path    = file.path(ROOT, "benchmarks", "stepstone_bench_data.parquet"),
    targets_path = file.path(ROOT, "benchmarks", "stepstone_bench_targets.json")
  ),
  small_hicard = list(
    data_path    = file.path(ROOT, "benchmarks", "yh0l", "small_hicard_data.parquet"),
    targets_path = file.path(ROOT, "benchmarks", "yh0l", "small_hicard_targets.json")
  )
)

# ── Metrics computation ───────────────────────────────────────────────────────
compute_metrics <- function(w, df, tgt) {
  W       <- sum(w)
  max_err <- 0.0
  mkl     <- 0.0
  chi2    <- 0.0
  for (v in names(tgt)) {
    agg <- tapply(w, df[[v]], sum)
    for (lv in names(tgt[[v]])) {
      S <- (agg[lv] %||% 0.0) / W
      T <- tgt[[v]][[lv]]
      max_err <- max(max_err, abs(S - T))
      if (T > 0 && S > 0) mkl <- mkl + T * log(T / S)
      chi2 <- chi2 + (S - T)^2 / max(T, 1e-12)
    }
  }
  wm  <- W / length(w)
  pos <- w > 0
  wkl <- sum(w[pos] * log(w[pos] / wm)) / W
  list(margin_kl = mkl, max_err = max_err, chi2 = chi2, weight_kl = wkl)
}

# ── Load a DGP ────────────────────────────────────────────────────────────────
load_dgp <- function(dgp) {
  df_raw  <- as.data.frame(read_parquet(dgp$data_path))
  tgt_raw <- fromJSON(dgp$targets_path)
  tgt     <- lapply(tgt_raw, function(x) { x <- unlist(x); x / sum(x) })
  df <- df_raw[, intersect(names(tgt), names(df_raw)), drop = FALSE]
  for (v in names(tgt)) {
    if (v %in% names(df)) df[[v]] <- factor(df[[v]], levels = names(tgt[[v]]))
  }
  list(df = df, tgt = tgt)
}

# ── Run one DGP ───────────────────────────────────────────────────────────────
run_dgp <- function(dgp_name, dgp) {
  cat(sprintf("\n=== R DGP: %s (chebyshev) ===\n", dgp_name))
  d <- load_dgp(dgp)
  df  <- d$df
  tgt <- d$tgt
  n   <- nrow(df)
  cat(sprintf("  n=%s margins=%d\n", format(n, big.mark = ","), length(tgt)))

  rows <- vector("list", NREPS)
  for (rep in seq_len(NREPS)) {
    cat(sprintf("  rep %d/%d ... ", rep, NREPS))
    t0 <- proc.time()[["elapsed"]]
    res <- harvest(
      df,
      tgt,
      method         = "chebyshev",
      convergence    = list(tol = TOL),
      max_weight     = MAX_WT,
      max_iterations = MAX_IT,
      attach_weights = FALSE,
      verbose        = 0L
    )
    t1      <- proc.time()[["elapsed"]]
    wall_ms <- (t1 - t0) * 1000.0
    cr      <- attr(res, "result")
    st      <- cr$status %||% -1L
    iters   <- cr$iterations %||% -1L
    me      <- cr$max_error %||% NA_real_
    cat(sprintf("status=%s iter=%d max_err=%.3e wall=%.0fms\n", st, iters, me, wall_ms))

    w_path <- file.path(WDIR, sprintf("%s_rep%d.feather", dgp_name, rep))
    write_feather(data.frame(w = as.double(res)), w_path)

    w <- as.double(res)
    m <- compute_metrics(w, df, tgt)

    rows[[rep]] <- data.frame(
      dgp        = dgp_name,
      rep        = rep,
      n          = n,
      margin_kl  = m$margin_kl,
      weight_kl  = m$weight_kl,
      max_err    = m$max_err,
      chi2       = m$chi2,
      wall_ms    = wall_ms,
      status     = st,
      iterations = iters,
      lang       = "R",
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

# ── Main ──────────────────────────────────────────────────────────────────────
all_results <- vector("list", length(dgps))
names(all_results) <- names(dgps)

for (dgp_name in names(dgps)) {
  dgp <- dgps[[dgp_name]]
  if (!file.exists(dgp$data_path)) {
    cat(sprintf("[SKIP] %s: data not found at %s\n", dgp_name, dgp$data_path))
    next
  }
  tryCatch({
    res_df <- run_dgp(dgp_name, dgp)
    out_path <- file.path(APMDIR, sprintf("results_R_%s.csv", dgp_name))
    write.csv(res_df, out_path, row.names = FALSE)
    cat(sprintf("  -> saved %s\n", out_path))
    all_results[[dgp_name]] <- res_df
  }, error = function(e) {
    cat(sprintf("[ERROR] DGP %s: %s\n", dgp_name, conditionMessage(e)))
  })
}

valid <- Filter(Negate(is.null), all_results)
if (length(valid) > 0) {
  agg <- do.call(rbind, valid)
  write.csv(agg, file.path(APMDIR, "results_R_all.csv"), row.names = FALSE)
  cat(sprintf("\nAll R results saved to benchmarks/2apm/results_R_all.csv\n"))
  print(agg[, c("dgp", "rep", "margin_kl", "weight_kl", "max_err", "wall_ms", "status", "iterations")])
}
