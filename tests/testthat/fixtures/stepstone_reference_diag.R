#!/usr/bin/env Rscript
# Diagnostic report on the captured stepstone_reference.rds.
# Computes per-margin fit metrics for BOTH solvers:
#   - max_err  = max_{k,j} |S_kj/W_total - target_kj|           (worst cell)
#   - L1_err   = Σ_{k,j} |S_kj/W_total - target_kj|              (summed per margin)
#   - chi2_kj  = (S_kj - target_kj·W_total)² / (target_kj·W_total)
#   - chi2     = Σ_{k,j} chi2_kj  (standard goodness-of-fit)
#   - KL_kj    = target_kj · log(target_kj / (S_kj/W_total))     (observed vs target)
#   - KL       = Σ_{k,j} KL_kj   (per-margin total; Kullback-Leibler divergence)
# Reports per-margin breakdown + aggregate; also extreme cells per solver.

suppressPackageStartupMessages({
  library(arrow); library(jsonlite)
})

FIXTURE_PATH <- "tests/testthat/fixtures/stepstone_reference.rds"
DATA_PATH    <- "benchmarks/stepstone_fulldata_bench_data.parquet"
TARGETS_PATH <- "benchmarks/stepstone_fulldata_bench_targets.json"

stopifnot(file.exists(FIXTURE_PATH), file.exists(DATA_PATH), file.exists(TARGETS_PATH))

ref <- readRDS(FIXTURE_PATH)
df  <- read_parquet(DATA_PATH)
if ("uuid" %in% names(df)) df$uuid <- NULL
tgt <- fromJSON(TARGETS_PATH)
tgt <- lapply(tgt, function(t) unlist(t))

# Apply same filtering as reference script (drop absent + zero-obs + renorm)
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
  per_margin <- lapply(names(tgt), function(nm) {
    lv <- names(tgt[[nm]])
    # Observed S_kj per category: sum of weights in each bucket
    S <- tapply(weights, df[[nm]], sum, default = 0)
    # Align to target order (some levels may be absent)
    S <- S[lv]; S[is.na(S)] <- 0
    target_kj <- unname(tgt[[nm]])
    S_rel     <- S / W_total              # observed proportion
    # Metrics per margin
    err_kj    <- abs(S_rel - target_kj)
    max_err   <- max(err_kj)
    L1        <- sum(err_kj)
    expect    <- target_kj * W_total
    chi2_kj   <- ifelse(expect > 0, (S - expect)^2 / expect, 0)
    chi2      <- sum(chi2_kj)
    # KL: target · log(target / observed); skip where target or observed is 0
    safe      <- target_kj > 0 & S_rel > 0
    KL_kj     <- ifelse(safe, target_kj * log(target_kj / pmax(S_rel, .Machine$double.eps)), 0)
    KL        <- sum(KL_kj)
    # Worst cell index
    worst_j   <- which.max(err_kj)
    list(margin = nm,
         n_cats = length(lv),
         max_err = max_err,
         L1_err  = L1,
         chi2    = chi2,
         KL      = KL,
         worst_cat = lv[worst_j],
         worst_err = err_kj[worst_j],
         worst_obs = S_rel[worst_j],
         worst_target = target_kj[worst_j])
  })
  # Aggregate
  max_err_all <- max(sapply(per_margin, `[[`, "max_err"))
  L1_all      <- sum(sapply(per_margin, `[[`, "L1_err"))
  chi2_all    <- sum(sapply(per_margin, `[[`, "chi2"))
  KL_all      <- sum(sapply(per_margin, `[[`, "KL"))
  list(per_margin = per_margin,
       max_err = max_err_all, L1_err = L1_all, chi2 = chi2_all, KL = KL_all,
       W_total = W_total, DEFF = length(weights) * sum(weights^2) / sum(weights)^2,
       ESS = sum(weights)^2 / sum(weights^2))
}

cat("=== Fit diagnostics: autumn vs leafblower ===\n")
cat(sprintf("n=%s K=%d total_cats=%d\n",
    format(length(ref$autumn_weights), big.mark=","), ref$K, ref$total_cats))

autumn_fit <- compute_fit(ref$autumn_weights, df, tgt)
lb_fit     <- compute_fit(ref$lb_weights,     df, tgt)

cat("\n=== Aggregate fit ===\n")
cat(sprintf("%-40s  %12s  %12s\n", "metric", "autumn", "leafblower"))
cat(sprintf("%-40s  %12s  %12s\n",
    strrep("-", 40), strrep("-", 12), strrep("-", 12)))
cat(sprintf("%-40s  %12.3e  %12.3e\n", "max_err (max over all cells)",
    autumn_fit$max_err, lb_fit$max_err))
cat(sprintf("%-40s  %12.3e  %12.3e\n", "L1_err (sum |S/W - target| all)",
    autumn_fit$L1_err, lb_fit$L1_err))
cat(sprintf("%-40s  %12.3e  %12.3e\n", "chi² (sum (S-E)²/E all margins)",
    autumn_fit$chi2, lb_fit$chi2))
cat(sprintf("%-40s  %12.3e  %12.3e\n", "KL divergence (target || observed)",
    autumn_fit$KL, lb_fit$KL))
cat(sprintf("%-40s  %12.3f  %12.3f\n", "DEFF", autumn_fit$DEFF, lb_fit$DEFF))
cat(sprintf("%-40s  %12s  %12s\n", "ESS",
    format(round(autumn_fit$ESS), big.mark=","),
    format(round(lb_fit$ESS),     big.mark=",")))

cat("\n=== Per-margin breakdown ===\n")
cat(sprintf("%-32s %4s | %-8s %-8s %-8s %-8s | %-8s %-8s %-8s %-8s\n",
    "margin", "cats",
    "aut_max", "aut_L1", "aut_chi2", "aut_KL",
    "lb_max",  "lb_L1",  "lb_chi2",  "lb_KL"))
cat(strrep("-", 150), "\n")
for (i in seq_along(ref$margin_order)) {
  nm <- ref$margin_order[i]
  a <- autumn_fit$per_margin[[i]]; l <- lb_fit$per_margin[[i]]
  cat(sprintf("%-32s %4d | %.2e %.2e %.2e %.2e | %.2e %.2e %.2e %.2e\n",
      nm, a$n_cats,
      a$max_err, a$L1_err, a$chi2, a$KL,
      l$max_err, l$L1_err, l$chi2, l$KL))
}

cat("\n=== Worst-cell per margin (larger max_err shown) ===\n")
cat(sprintf("%-32s %-40s  %-11s %-10s %-10s\n",
    "margin", "worst_cat (autumn | leafblower)", "err(aut/lb)", "obs", "target"))
for (i in seq_along(ref$margin_order)) {
  nm <- ref$margin_order[i]
  a <- autumn_fit$per_margin[[i]]; l <- lb_fit$per_margin[[i]]
  larger <- if (a$max_err > l$max_err) "autumn" else "leafblower"
  worst <- if (larger == "autumn") a else l
  cat(sprintf("  %-30s %-40s  %.2e  %.4f  %.4f  [%s]\n",
      nm, substr(worst$worst_cat, 1, 40), worst$worst_err,
      worst$worst_obs, worst$worst_target, larger))
}

cat("\n=== Relative comparison ===\n")
cat(sprintf("  autumn / leafblower  max_err: %.2fx  L1: %.2fx  chi²: %.2fx  KL: %.2fx\n",
    autumn_fit$max_err / lb_fit$max_err,
    autumn_fit$L1_err  / lb_fit$L1_err,
    autumn_fit$chi2    / lb_fit$chi2,
    autumn_fit$KL      / lb_fit$KL))
