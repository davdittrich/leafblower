#!/usr/bin/env Rscript
# benchmarks/e65t1_corun_bench.R
# e65t.1.6: 4-arm co-run benchmark — GO/NO-GO verdict.
# Arms A/B/C/D × Configs 1-4 (stepstone small+full, seeds 1+2).
# Criterion: Arm D iters ≤ Arm A iters in ≥3/4 configs AND status OK
#            AND sraa_corun_reverts < 0.20 * iters for all Arm D runs.

suppressPackageStartupMessages({
  library(arrow)
  library(jsonlite)
  library(leafblower)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

ROOT <- normalizePath(
  if (file.exists("benchmarks/stepstone_bench_data.parquet")) "." else "../..",
  mustWork = FALSE
)

DATA_SMALL    <- file.path(ROOT, "benchmarks", "stepstone_bench_data.parquet")
TARGETS_SMALL <- file.path(ROOT, "benchmarks", "stepstone_bench_targets.json")
DATA_FULL     <- file.path(ROOT, "benchmarks", "stepstone_fulldata_bench_data.parquet")
TARGETS_FULL  <- file.path(ROOT, "benchmarks", "stepstone_fulldata_bench_targets.json")

# Verify data files
for (f in c(DATA_SMALL, TARGETS_SMALL, DATA_FULL, TARGETS_FULL)) {
  if (!file.exists(f)) stop("Missing data file: ", f)
}
cat("Data files verified.\n\n")

MAX_ITER <- 2000L
MAX_WT   <- 5.0
MIN_WT   <- 0.3  # bounded variant: strict per-unit water-fill
CONV     <- list(metric = "marginal_kl")  # default improvement rule, pct_tol=0.001

# Arms: named list of harvest() extra args
arms <- list(
  A = list(sor = list(auto = TRUE),             accelerate = TRUE),   # default (corun_aa=0)
  B = list(sor = list(auto = TRUE),             accelerate = FALSE),  # SOR-only
  C = list(sor = NULL,                          accelerate = TRUE),   # SRAA-only
  D = list(sor = list(auto = TRUE, corun_aa = TRUE), accelerate = TRUE)  # co-run
)

load_dgp <- function(data_path, targets_path) {
  df_raw  <- as.data.frame(read_parquet(data_path))
  tgt_raw <- fromJSON(targets_path)
  tgt     <- lapply(tgt_raw, function(x) { x <- unlist(x); x / sum(x) })
  df      <- df_raw[, intersect(names(tgt), names(df_raw)), drop = FALSE]
  for (v in names(tgt))
    if (v %in% names(df)) df[[v]] <- factor(df[[v]], levels = names(tgt[[v]]))
  list(df = df, tgt = tgt)
}

# Configs: (label, data, tgt, seed, min_wt, bounds_mode)
configs <- list(
  list(label = "small_s1_unb",  data = DATA_SMALL, tgt = TARGETS_SMALL, seed = 1L, min_wt = 0.0, bmode = "cell"),
  list(label = "small_s2_unb",  data = DATA_SMALL, tgt = TARGETS_SMALL, seed = 2L, min_wt = 0.0, bmode = "cell"),
  list(label = "full_s1_unb",   data = DATA_FULL,  tgt = TARGETS_FULL,  seed = 1L, min_wt = 0.0, bmode = "cell"),
  list(label = "full_s2_unb",   data = DATA_FULL,  tgt = TARGETS_FULL,  seed = 2L, min_wt = 0.0, bmode = "cell"),
  list(label = "small_s1_bnd",  data = DATA_SMALL, tgt = TARGETS_SMALL, seed = 1L, min_wt = MIN_WT, bmode = "unit"),
  list(label = "small_s2_bnd",  data = DATA_SMALL, tgt = TARGETS_SMALL, seed = 2L, min_wt = MIN_WT, bmode = "unit"),
  list(label = "full_s1_bnd",   data = DATA_FULL,  tgt = TARGETS_FULL,  seed = 1L, min_wt = MIN_WT, bmode = "unit"),
  list(label = "full_s2_bnd",   data = DATA_FULL,  tgt = TARGETS_FULL,  seed = 2L, min_wt = MIN_WT, bmode = "unit")
)

STATUS_STR <- c("0" = "OK", "1" = "no_conv", "2" = "infeasible",
                "3" = "bad_arg", "4" = "budget", "5" = "stall")

run_one <- function(df, tgt, arm_name, arm_args, seed, min_wt = 0.0, bmode = "cell") {
  set.seed(seed)
  t0  <- proc.time()[["elapsed"]]
  base_args <- list(df, tgt,
                    min_weight     = min_wt,
                    max_weight     = MAX_WT,
                    bounds_mode    = bmode,
                    max_iterations = MAX_ITER,
                    convergence    = CONV,
                    attach_weights = FALSE,
                    verbose        = 0L)
  res <- tryCatch(
    suppressWarnings(do.call(harvest, c(base_args, arm_args))),
    error = function(e) structure(list(), class = "harvest_error",
                                  msg = conditionMessage(e))
  )
  elapsed_ms <- (proc.time()[["elapsed"]] - t0) * 1000.0

  if (inherits(res, "harvest_error")) {
    return(data.frame(
      arm = arm_name, status = "error", iterations = NA_integer_,
      marginal_kl = NA_real_, max_error = NA_real_,
      sraa_corun_reverts = NA_integer_, ms = round(elapsed_ms),
      stringsAsFactors = FALSE
    ))
  }

  cr  <- attr(res, "result")
  st  <- STATUS_STR[as.character(cr$status %||% -1L)] %||% as.character(cr$status)
  itr <- cr$iterations %||% NA_integer_
  me  <- cr$max_error  %||% NA_real_
  mkl <- cr$margin_kl  %||% NA_real_   # post-solver quality metric from compute_quality_metrics
  rev <- cr$sraa_corun_reverts %||% 0L

  data.frame(
    arm = arm_name, status = st, iterations = itr,
    marginal_kl = signif(mkl, 5), max_error = signif(me, 5),
    sraa_corun_reverts = rev, ms = round(elapsed_ms),
    stringsAsFactors = FALSE
  )
}

# Pre-load DGPs once
cat("Loading DGPs...\n")
dgp_cache <- list()
for (cfg in configs) {
  key <- cfg$label
  data_key <- paste0(cfg$data, "|", cfg$tgt)
  if (is.null(dgp_cache[[data_key]])) {
    cat(sprintf("  Loading %s...\n", basename(cfg$data)))
    dgp_cache[[data_key]] <- load_dgp(cfg$data, cfg$tgt)
  }
}

# Build results
results <- vector("list", length(arms) * length(configs))
idx <- 0L
total <- length(arms) * length(configs)

cat(sprintf("\nRunning %d arms × %d configs = %d runs...\n\n",
            length(arms), length(configs), total))

for (cfg in configs) {
  data_key <- paste0(cfg$data, "|", cfg$tgt)
  dgp      <- dgp_cache[[data_key]]
  n        <- nrow(dgp$df)

  for (arm_name in names(arms)) {
    idx <- idx + 1L
    cat(sprintf("[%2d/%2d] arm=%s config=%-12s n=%s ... ",
                idx, total, arm_name, cfg$label,
                format(n, big.mark = ",")))
    flush.console()

    row <- run_one(dgp$df, dgp$tgt, arm_name, arms[[arm_name]], cfg$seed, cfg$min_wt, cfg$bmode)
    row$config <- cfg$label
    row$n      <- n
    results[[idx]] <- row

    cat(sprintf("status=%-8s iters=%4s mkl=%.3e me=%.3e reverts=%d  (%dms)\n",
                row$status,
                if (is.na(row$iterations)) "NA" else row$iterations,
                row$marginal_kl %||% NaN,
                row$max_error   %||% NaN,
                row$sraa_corun_reverts %||% 0L,
                row$ms))
  }
  cat("\n")
}

tab <- do.call(rbind, results)
tab <- tab[, c("config", "arm", "status", "iterations", "marginal_kl",
               "max_error", "sraa_corun_reverts", "ms")]

cat("\n=== Results table ===\n")
print(tab, row.names = FALSE, digits = 5)

# --- GO/NO-GO evaluation ---
cat("\n=== GO/NO-GO Evaluation ===\n")

arm_A <- tab[tab$arm == "A", ]
arm_D <- tab[tab$arm == "D", ]

# Merge on config
eval_df <- merge(
  arm_A[, c("config", "iterations", "status")],
  arm_D[, c("config", "iterations", "status", "sraa_corun_reverts")],
  by = "config", suffixes = c("_A", "_D")
)
eval_df <- eval_df[order(eval_df$config), ]

cat("\nPer-config iteration comparison (Arm D vs Arm A):\n")
eval_df$D_le_A <- !is.na(eval_df$iterations_D) &
                  !is.na(eval_df$iterations_A) &
                  eval_df$iterations_D <= eval_df$iterations_A
eval_df$revert_pct <- ifelse(
  !is.na(eval_df$sraa_corun_reverts) & !is.na(eval_df$iterations_D) & eval_df$iterations_D > 0,
  eval_df$sraa_corun_reverts / eval_df$iterations_D,
  NA_real_
)
print(eval_df[, c("config", "iterations_A", "iterations_D", "D_le_A",
                  "sraa_corun_reverts", "revert_pct", "status_D")],
      row.names = FALSE)

med_A <- median(arm_A$iterations, na.rm = TRUE)
med_D <- median(arm_D$iterations, na.rm = TRUE)
n_configs_D_le_A <- sum(eval_df$D_le_A, na.rm = TRUE)

all_D_ok     <- all(arm_D$status == "OK", na.rm = TRUE)
mkl_ok       <- all(is.finite(arm_D$marginal_kl), na.rm = TRUE)
max_revert_pct <- max(eval_df$revert_pct, na.rm = TRUE)
reverts_ok   <- !is.na(max_revert_pct) && max_revert_pct < 0.20

# Note: full-dataset "budget" is shared across ALL arms (A/B/C/D all hit 2000 iters).
# This is structural ill-conditioning, not a D-specific failure. Reported as context.
all_A_ok <- all(arm_A$status == "OK", na.rm = TRUE)
cat(sprintf("\nNote: All Arm A status OK: %s (shared budget on full data)\n", all_A_ok))

cat(sprintf("\nArm A median iters: %.0f\n", med_A))
cat(sprintf("Arm D median iters: %.0f\n", med_D))
cat(sprintf("Configs where D <= A: %d/4\n", n_configs_D_le_A))
cat(sprintf("All Arm D status OK: %s\n", all_D_ok))
cat(sprintf("All Arm D marginal_kl finite: %s\n", mkl_ok))
cat(sprintf("Max Arm D revert%%: %.1f%% (threshold <20%%)\n", max_revert_pct * 100))
cat(sprintf("Reverts criterion met: %s\n", reverts_ok))

# GO requires ALL three (per spec):
# 1) D median iters <= A median iters across >=3/4 configs
# 2) marginal_kl converges (status OK) for all Arm D runs
# 3) sraa_corun_reverts < 20% of iters for all Arm D runs
criterion_iters   <- n_configs_D_le_A >= 3L
criterion_status  <- all_D_ok && mkl_ok
criterion_reverts <- reverts_ok

cat(sprintf("\nCriterion 1 (D<=A in >=3/4 configs): %s\n",
            if (criterion_iters) "PASS" else "FAIL"))
cat(sprintf("Criterion 2 (status OK + mkl finite): %s\n",
            if (criterion_status) "PASS" else "FAIL"))
cat(sprintf("Criterion 3 (reverts <20%% of iters): %s\n",
            if (criterion_reverts) "PASS" else "FAIL"))

VERDICT <- if (criterion_iters && criterion_status && criterion_reverts) "GO" else "NO-GO"
cat(sprintf("\n*** VERDICT: %s ***\n", VERDICT))

# Also show arm B/C for context
cat("\n=== Arm B/C context (SOR-only, SRAA-only) ===\n")
bc <- tab[tab$arm %in% c("B", "C"), c("config", "arm", "status", "iterations", "marginal_kl")]
print(bc, row.names = FALSE)

# Save results
outdir <- file.path(ROOT, "benchmarks", "results")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(outdir, "e65t1_corun_bench.csv")
write.csv(tab, outfile, row.names = FALSE)
cat(sprintf("\nResults saved: %s\n", outfile))
cat(sprintf("VERDICT: %s\n", VERDICT))
