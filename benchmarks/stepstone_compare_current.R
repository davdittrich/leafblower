#!/usr/bin/env Rscript
# 3-way stepstone comparison: autumn (ref) vs leafblower cell-mode (ref+current)
# vs leafblower unit-mode (current, P3.1 strict per-obs bounds).
# autumn unchanged across HEAD revs; leafblower rerun at current HEAD.

Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(arrow); library(jsonlite); library(leafblower)
})

FIXTURE  <- "tests/testthat/fixtures/stepstone_reference.rds"
DATA_P   <- "benchmarks/stepstone_fulldata_bench_data.parquet"
TARGETS_P<- "benchmarks/stepstone_fulldata_bench_targets.json"

ref <- readRDS(FIXTURE)
df  <- read_parquet(DATA_P)
if ("uuid" %in% names(df)) df$uuid <- NULL
tgt <- fromJSON(TARGETS_P)
tgt <- lapply(tgt, function(t) unlist(t))
for (nm in names(tgt)) {
  sv <- table(df[[nm]]); present <- names(sv)[sv > 0]
  t_ <- tgt[[nm]]; keep <- names(t_) %in% present; t_ <- t_[keep]
  n_per_cat <- as.integer(sv[names(t_)])
  zero_obs <- which(n_per_cat == 0 & t_ > 0)
  if (length(zero_obs) > 0) t_ <- t_[-zero_obs]
  tgt[[nm]] <- t_ / sum(t_)
}

compute_fit <- function(weights, df, tgt) {
  W_total <- sum(weights)
  per_margin_err <- c(); per_margin_L1 <- 0
  per_margin_chi2 <- 0;  per_margin_KL <- 0
  for (nm in names(tgt)) {
    lv <- names(tgt[[nm]])
    S  <- tapply(weights, df[[nm]], sum, default=0)[lv]
    S[is.na(S)] <- 0
    tk <- unname(tgt[[nm]])
    Sr <- S / W_total
    err_kj <- abs(Sr - tk)
    per_margin_err  <- c(per_margin_err, max(err_kj))
    per_margin_L1   <- per_margin_L1 + sum(err_kj)
    expect <- tk * W_total
    per_margin_chi2 <- per_margin_chi2 + sum(ifelse(expect>0, (S-expect)^2/expect, 0))
    safe <- tk > 0 & Sr > 0
    per_margin_KL <- per_margin_KL + sum(ifelse(safe, tk*log(tk/pmax(Sr, .Machine$double.eps)), 0))
  }
  list(max_err = max(per_margin_err),
       L1      = per_margin_L1,
       chi2    = per_margin_chi2,
       KL      = per_margin_KL,
       DEFF    = length(weights)*sum(weights^2)/sum(weights)^2,
       ESS     = sum(weights)^2/sum(weights^2),
       wmin    = min(weights), wmed = median(weights), wmax = max(weights),
       n_over_5 = sum(weights > 5.0 + 1e-9),
       n_under_0 = sum(weights < 0 - 1e-9))
}

cat(sprintf("=== Reference SHA: %s ===\n", ref$git_sha))
cat(sprintf("=== Current HEAD:   %s ===\n",
    system("git rev-parse HEAD", intern=TRUE)))

# Cell-mode (default)
cat("\n=== Current lb ORIS — cell-mode (bounds_mode='cell', default) ===\n")
t0 <- Sys.time()
r_cell <- suppressWarnings(
  leafblower::harvest(df, tgt, method="oris",
    max_weight=5, min_weight=0,
    max_iterations=3000L, convergence=list(absolute=1e-10),
    bounds_mode="cell", attach_weights=FALSE)
)
wall_cell <- as.numeric(Sys.time()-t0, units="secs")
w_cell  <- as.numeric(r_cell)
info_cell <- attr(r_cell, "result")
cat(sprintf("  wall=%.1fs status=%d iter=%d solver_err=%.3e n_bounds_violated=%d n_bounds_clamped=%d\n",
    wall_cell, info_cell$status, info_cell$iterations, info_cell$max_error,
    info_cell$n_bounds_violated, info_cell$n_bounds_clamped))

# Unit-mode (P3.1 strict per-obs bounds)
cat("\n=== Current lb ORIS — unit-mode (bounds_mode='unit', P3.1 strict) ===\n")
t0 <- Sys.time()
r_unit <- suppressWarnings(
  leafblower::harvest(df, tgt, method="oris",
    max_weight=5, min_weight=0,
    max_iterations=3000L, convergence=list(absolute=1e-10),
    bounds_mode="unit", attach_weights=FALSE)
)
wall_unit <- as.numeric(Sys.time()-t0, units="secs")
w_unit  <- as.numeric(r_unit)
info_unit <- attr(r_unit, "result")
cat(sprintf("  wall=%.1fs status=%d iter=%d solver_err=%.3e n_bounds_violated=%d n_bounds_clamped=%d\n",
    wall_unit, info_unit$status, info_unit$iterations, info_unit$max_error,
    info_unit$n_bounds_violated, info_unit$n_bounds_clamped))

fit_autumn  <- compute_fit(ref$autumn_weights, df, tgt)
fit_lb_ref  <- compute_fit(ref$lb_weights,     df, tgt)
fit_cell    <- compute_fit(w_cell,             df, tgt)
fit_unit    <- compute_fit(w_unit,             df, tgt)

cat("\n=== 4-way fit comparison ===\n")
cat(sprintf("%-15s  %14s  %14s  %14s  %14s\n", "metric",
    "autumn(ref)", "lb cell(ref)", "lb cell(cur)", "lb unit(cur)"))
cat(sprintf("%-15s  %14s  %14s  %14s  %14s\n",
    strrep("-",15), strrep("-",14), strrep("-",14), strrep("-",14), strrep("-",14)))

fmt <- function(label, a, b, c, d, f="%14.3e") {
  cat(sprintf(paste0("%-15s  ", f, "  ", f, "  ", f, "  ", f, "\n"), label, a, b, c, d))
}
fmt("max_err",  fit_autumn$max_err, fit_lb_ref$max_err, fit_cell$max_err, fit_unit$max_err)
fmt("L1_err",   fit_autumn$L1,      fit_lb_ref$L1,      fit_cell$L1,      fit_unit$L1)
fmt("chi²",    fit_autumn$chi2,    fit_lb_ref$chi2,    fit_cell$chi2,    fit_unit$chi2)
fmt("KL_div",   fit_autumn$KL,      fit_lb_ref$KL,      fit_cell$KL,      fit_unit$KL)
fmt("DEFF",     fit_autumn$DEFF,    fit_lb_ref$DEFF,    fit_cell$DEFF,    fit_unit$DEFF, "%14.4f")
cat(sprintf("%-15s  %14s  %14s  %14s  %14s\n", "ESS",
    format(round(fit_autumn$ESS), big.mark=","),
    format(round(fit_lb_ref$ESS), big.mark=","),
    format(round(fit_cell$ESS),   big.mark=","),
    format(round(fit_unit$ESS),   big.mark=",")))
fmt("w_min",    fit_autumn$wmin, fit_lb_ref$wmin, fit_cell$wmin, fit_unit$wmin, "%14.4f")
fmt("w_median", fit_autumn$wmed, fit_lb_ref$wmed, fit_cell$wmed, fit_unit$wmed, "%14.4f")
fmt("w_max",    fit_autumn$wmax, fit_lb_ref$wmax, fit_cell$wmax, fit_unit$wmax, "%14.4f")
cat(sprintf("%-15s  %14d  %14d  %14d  %14d\n", "w > 5.0",
    fit_autumn$n_over_5, fit_lb_ref$n_over_5, fit_cell$n_over_5, fit_unit$n_over_5))
cat(sprintf("%-15s  %14d  %14d  %14d  %14d\n", "w < 0.0",
    fit_autumn$n_under_0, fit_lb_ref$n_under_0, fit_cell$n_under_0, fit_unit$n_under_0))

cat("\n=== P3.1 unit-mode strict-bounds acceptance ===\n")
cat(sprintf("  max(w) ≤ 5.0+1e-9:   %s (observed max = %.10f)\n",
    fit_unit$wmax <= 5.0 + 1e-9, fit_unit$wmax))
cat(sprintf("  min(w) ≥ 0.0-1e-9:   %s (observed min = %.10f)\n",
    fit_unit$wmin >= -1e-9, fit_unit$wmin))
cat(sprintf("  n_bounds_clamped:    %d (post-leafblower-kssd: running counter reports every clamp)\n",
    info_unit$n_bounds_clamped))

cat("\n=== Regression: cell-mode current vs leafblower reference ===\n")
dw <- abs(w_cell - ref$lb_weights)
cat(sprintf("  max |Δw|: %.3e   Pearson r: %.6f\n",
    max(dw), cor(w_cell, ref$lb_weights)))

cat("\n=== Correlations ===\n")
cat(sprintf("  autumn      ↔ lb cell(cur): r=%.6f\n", cor(ref$autumn_weights, w_cell)))
cat(sprintf("  autumn      ↔ lb unit(cur): r=%.6f\n", cor(ref$autumn_weights, w_unit)))
cat(sprintf("  lb cell(cur)↔ lb unit(cur): r=%.6f\n", cor(w_cell, w_unit)))
