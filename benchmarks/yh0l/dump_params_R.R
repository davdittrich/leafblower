#!/usr/bin/env Rscript
# T1: dump all 37 C_rk_calibrate args for ieppa_soft as R resolves them.
# Call: harvest(method='ieppa_soft', convergence=list(tol=1e-4),
#               max_weight=5, max_iterations=3000)
# Does NOT call the solver; replicates harvest.R's arg-prep then prints.

if (!file.exists("benchmarks/yh0l")) {
  stop("Run from project root. Current wd=", getwd())
}

suppressPackageStartupMessages({
  library(leafblower)
})

# Replicate harvest.R logic for:
#   method         = "ieppa_soft"
#   convergence    = list(tol = 1e-4)
#   max_weight     = 5
#   max_iterations = 3000
#   all other args = their defaults

method          <- "ieppa_soft"
convergence_arg <- list(tol = 1e-4)
min_weight      <- 0.0
max_weight      <- 5.0
max_iterations  <- 3000L
verbose         <- 0L
capacity_penalty <- NULL
alm_penalty     <- NULL
bounds_mode     <- "cell"
homotopy_levels       <- 1L
homotopy_start_factor <- 1.0
homotopy_end_factor   <- 1.0
homotopy_budget_p     <- 0.5
scheduler       <- "round_robin"
eta_schedule    <- "fixed"
eta_start       <- 1.0
eta_end         <- 1.0
eta_schedule_power <- 0.5
accelerate      <- FALSE
newton_tsvd_ratio <- 1e-8
ridge_lambda    <- 0.0
sor             <- NULL   # R default: sor=NULL → disabled (parse_sor returns enabled=0)

# Replicate parse_convergence (harvest.R:758)
# convergence=list(tol=1e-4) → tol overridden, metric=default→marginal_kl(ieppa_soft), rule=improvement
conv_env <- new.env(parent = emptyenv())
environment(leafblower:::parse_convergence) -> orig_env
conv <- leafblower:::parse_convergence(convergence_arg)

# Per-method metric override (harvest.R:308-331): no explicit metric → ieppa_soft = marginal_kl
conv$metric <- "marginal_kl"

# parse_sor: sor=NULL → disabled
sor_cfg <- leafblower:::parse_sor(sor)

# bounds_mode
bounds_mode_int <- match("cell", c("cell", "unit")) - 1L  # 0

# metric/rule/stop_when int mapping
metric_int    <- c(max_err=0L, mean_err=1L, kl=2L, chi2=3L,
                   grake_norm=4L, l1_weight=5L, marginal_kl=6L)
rule_int      <- c(threshold=0L, improvement=1L, plateau=2L)
stop_when_int <- c(any=0L, all=1L)

# tol_abs (legacy): absolute_tol > 0 → use it; else 1e-6
# conv$absolute_tol comes from parse_convergence; for list(tol=1e-4) this is 0.0
# (rule="improvement" → pct_tol=tol=1e-4, absolute_tol=0.0)
tol_abs_val <- if (conv$absolute_tol > 0) conv$absolute_tol else 1e-6

# scheduler/eta_schedule → character (passed as string to C bridge)
scheduler_str    <- "round_robin"
eta_schedule_str <- "fixed"

# accelerate_bool
accelerate_bool <- FALSE

cat("=== R arg dump: harvest(method='ieppa_soft', convergence=list(tol=1e-4), max_weight=5, max_iterations=3000) ===\n\n")

# Slots 1-3: data-dependent structures (structural summary only)
cat("  slot  1 | group_ids_r            = VECSXP[K]; each INTSXP len=n (0-idx codes, -1=NA). K=n_margins, n=nrow(data)\n")
cat("  slot  2 | cat_counts_r           = INTSXP[K]; per-margin level counts. sum=total_cells\n")
cat("  slot  3 | targets_r              = VECSXP[K]; each REALSXP len=cat_counts[k]. sum(targets[[k]])≈1\n")

args_list <- list(
  n_obs                = "n_obs (INTSXP scalar = nrow(data))",
  min_weight           = as.double(min_weight),
  max_weight           = as.double(max_weight),
  method               = as.character(method),
  verbose              = as.integer(verbose),
  max_iterations       = as.integer(max_iterations),
  start_weights        = "sw_vec (REALSXP len=n or NULL→uniform 1.0)",
  # slot 11
  capacity_penalty     = if (is.null(capacity_penalty)) -1.0 else as.double(capacity_penalty),
  # slot 10
  alm_penalty          = if (is.null(alm_penalty)) -1.0 else as.double(alm_penalty),
  # slot 11: legacy tol_abs
  tol_abs              = as.double(tol_abs_val),
  bounds_mode          = as.integer(bounds_mode_int),
  homotopy_levels      = as.integer(homotopy_levels),
  homotopy_start_factor = as.double(homotopy_start_factor),
  homotopy_end_factor  = as.double(homotopy_end_factor),
  homotopy_budget_p    = as.double(homotopy_budget_p),
  scheduler            = as.character(scheduler_str),
  eta_schedule         = as.character(eta_schedule_str),
  eta_start            = as.double(eta_start),
  eta_end              = as.double(eta_end),
  eta_schedule_power   = as.double(eta_schedule_power),
  # WU-A convergence
  pct_tol              = as.double(conv$pct_tol),
  absolute_tol         = as.double(conv$absolute_tol),
  metric               = as.integer(metric_int[[conv$metric]]),
  rule                 = as.integer(rule_int[[conv$rule]]),
  stop_when            = as.integer(stop_when_int[[conv$stop_when]]),
  # SOR config
  sor_enabled          = as.integer(sor_cfg$enabled),
  sor_auto             = as.integer(sor_cfg$auto),
  sor_omega_init       = as.double(sor_cfg$omega_init),
  sor_omega_min        = as.double(sor_cfg$omega_min),
  sor_omega_fixed      = as.double(sor_cfg$omega_fixed),
  sor_burnin           = as.integer(sor_cfg$burnin),
  # SRAA
  accelerate           = as.integer(accelerate_bool),
  newton_tsvd_ratio    = as.double(newton_tsvd_ratio),
  ridge_lambda         = as.double(ridge_lambda)
)

# Print slots 4-37 (data/tuning args); slots 1-3 printed above as structural summaries
scalar_slot <- 4L  # slot numbering from harvest.R .Call
for (nm in names(args_list)) {
  cat(sprintf("  slot %2d | %-22s = %s\n", scalar_slot, nm, deparse(args_list[[nm]])))
  scalar_slot <- scalar_slot + 1L
}

cat("\n--- conv struct ---\n")
cat(sprintf("  conv$metric       = %s\n", conv$metric))
cat(sprintf("  conv$rule         = %s\n", conv$rule))
cat(sprintf("  conv$pct_tol      = %s\n", conv$pct_tol))
cat(sprintf("  conv$absolute_tol = %s\n", conv$absolute_tol))
cat(sprintf("  conv$stop_when    = %s\n", conv$stop_when))

cat("\n--- sor_cfg struct ---\n")
cat(sprintf("  enabled   = %d\n", sor_cfg$enabled))
cat(sprintf("  auto      = %d\n", sor_cfg$auto))
cat(sprintf("  omega_init= %.4f\n", sor_cfg$omega_init))
cat(sprintf("  omega_min = %.4f\n", sor_cfg$omega_min))
cat(sprintf("  omega_fixed= %.4f\n", sor_cfg$omega_fixed))
cat(sprintf("  burnin    = %d\n", sor_cfg$burnin))

cat("\nDone.\n")
