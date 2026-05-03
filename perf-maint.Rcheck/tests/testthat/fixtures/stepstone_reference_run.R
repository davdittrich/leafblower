#!/usr/bin/env Rscript
# Reference stepstone run matching the ORIGINAL Rmd
# (Salary-Data/code/22-weighting-create-weights-2025-firmsize.Rmd lines 468-476):
#   harvest(data = stepstone_imputed, target = target_anes,
#           max_iterations = 3000, max_weight = 5,
#           select_function = "all", error_function = "mean",
#           convergence = c(pct=1e-15, absolute=1e-10, time=NULL, single_weight=NULL))
#
# Captures margin counts, autumn and leafblower weights, DEFF/ESS for
# future benchmark comparison. Loads from the already-saved parquet;
# does NOT rebuild from /home/dd/stepstone/.
#
# Rmd expected (lines 484, 494 comments):
#   n = 1,248,521  ESS = 760,219.5  →  DEFF ≈ 1.64

Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({
  library(arrow); library(jsonlite); library(leafblower); library(autumn)
})

DATA_PATH    <- "benchmarks/stepstone_fulldata_bench_data.parquet"
TARGETS_PATH <- "benchmarks/stepstone_fulldata_bench_targets.json"
OUT_DIR      <- "tests/testthat/fixtures"
dir.create(OUT_DIR, recursive=TRUE, showWarnings=FALSE)

stopifnot(file.exists(DATA_PATH), file.exists(TARGETS_PATH))

df  <- read_parquet(DATA_PATH)
tgt <- fromJSON(TARGETS_PATH)
tgt <- lapply(tgt, function(t) unlist(t))

# Verify margin order matches original Rmd line 457-466
expected_order <- c(
  "rk_i_loc_time_age10_gender", "rk_i_loc_wz", "rk_age10_gender",
  "rk_i_loc_time_gender", "rk_i_loc_fsize", "rk_i_loc_Abschluss_gender",
  "rk_gender_time", "rk_time", "rk_gender"
)
cat("=== Margin verification vs Rmd line 456-466 ===\n")
if (!identical(names(tgt), expected_order)) {
  cat("MISMATCH:\n  got:      ", paste(names(tgt), collapse=", "), "\n")
  cat("  expected: ", paste(expected_order, collapse=", "), "\n")
  stop("margin order diverges from original Rmd")
}
cat("  OK (9 margins, correct order)\n")

if ("uuid" %in% names(df)) df$uuid <- NULL

# Drop target cells absent in survey + renormalize per margin
for (nm in names(tgt)) {
  sv <- table(df[[nm]])
  present <- names(sv)[sv > 0]
  t_ <- tgt[[nm]]
  keep <- names(t_) %in% present
  if (sum(!keep) > 0) t_ <- t_[keep]
  n_per_cat <- as.integer(sv[names(t_)])
  zero_obs <- which(n_per_cat == 0 & t_ > 0)
  if (length(zero_obs) > 0) t_ <- t_[-zero_obs]
  tgt[[nm]] <- t_ / sum(t_)
}

margin_counts <- sapply(tgt, length)
cat(sprintf("\n=== Final margin sizes ===\nn=%s K=%d total_cats=%d\n",
    format(nrow(df), big.mark=","), length(tgt), sum(margin_counts)))
for (nm in names(tgt)) cat(sprintf("  %s: %d cats\n", nm, length(tgt[[nm]])))

# Original Rmd parameters
MAX_ITER    <- 3000L
MAX_WEIGHT  <- 5
AUTUMN_CONV <- c(pct=1e-15, absolute=1e-10, time=NULL, single_weight=NULL)

# =============================================================
# 1. autumn::harvest — the ORIGINAL reference (matches Rmd exactly)
# =============================================================
cat("\n=== autumn::harvest (ORIGINAL Rmd params, max_iter=3000 abs=1e-10, accelerate=TRUE) ===\n")
t0 <- Sys.time()
r_autumn <- suppressWarnings(
  autumn::harvest(df, tgt,
                  max_iterations = MAX_ITER,
                  max_weight = MAX_WEIGHT,
                  select_function = "all",
                  error_function = "mean",
                  convergence = AUTUMN_CONV,
                  accelerate = TRUE)
)
autumn_wall_s <- as.numeric(Sys.time() - t0, units="secs")
w_autumn <- r_autumn$weights
cat(sprintf("  wall:      %.1fs\n", autumn_wall_s))
cat(sprintf("  n:         %s\n", format(length(w_autumn), big.mark=",")))
cat(sprintf("  w:         min=%.4f med=%.4f max=%.4f\n",
    min(w_autumn), median(w_autumn), max(w_autumn)))
cat(sprintf("  w > 5:     %d (%.2f%%)  [autumn: no strict box]\n",
    sum(w_autumn > 5), 100*mean(w_autumn > 5)))
autumn_deff <- autumn::design_effect(r_autumn)
autumn_ess  <- autumn::effective_sample_size(r_autumn)
cat(sprintf("  DEFF:      %.3f\n", autumn_deff))
cat(sprintf("  ESS:       %s  (Rmd comment: 760,219.5 → delta = %.1f)\n",
    format(round(autumn_ess), big.mark=","),
    autumn_ess - 760219.5))

# Save autumn intermediate immediately so leafblower failure doesn't lose it
intermediate <- list(
  autumn_wall_s = autumn_wall_s, autumn_weights = w_autumn,
  autumn_deff = autumn_deff, autumn_ess = autumn_ess
)
saveRDS(intermediate, file.path(OUT_DIR, "stepstone_reference_autumn_only.rds"))
cat("  autumn intermediate saved.\n")

# =============================================================
# 2. leafblower::harvest (ieppa) — same params, for comparison
# Qualified explicitly — autumn also exports harvest().
# =============================================================
cat("\n=== leafblower::harvest (ieppa, same params for comparison) ===\n")
t0 <- Sys.time()
r_lb <- suppressWarnings(
  leafblower::harvest(df, tgt, method="ieppa",
          max_weight=MAX_WEIGHT, min_weight=0,
          max_iterations=MAX_ITER,
          convergence=list(absolute=1e-10),
          attach_weights=FALSE)
)
lb_wall_s <- as.numeric(Sys.time() - t0, units="secs")
w_lb <- as.numeric(r_lb)
info_lb <- attr(r_lb, "result")
cat(sprintf("  wall:       %.1fs\n", lb_wall_s))
cat(sprintf("  status:     %d (0=OK, 1=NOCONV, 2=INFEAS)\n", info_lb$status))
cat(sprintf("  iterations: %d\n", info_lb$iterations))
cat(sprintf("  max_error:  %.3e\n", info_lb$max_error))
cat(sprintf("  w:          min=%.4f med=%.4f max=%.4f\n",
    min(w_lb), median(w_lb), max(w_lb)))
cat(sprintf("  w > 5:      %d (%.2f%%)  [leafblower: strict Dykstra box]\n",
    sum(w_lb > 5), 100*mean(w_lb > 5)))
lb_deff <- leafblower::design_effect(w_lb)
lb_ess  <- leafblower::effective_sample_size(w_lb)
cat(sprintf("  DEFF:       %.3f\n", lb_deff))
cat(sprintf("  ESS:        %s\n", format(round(lb_ess), big.mark=",")))

# =============================================================
# 3. Weight agreement autumn <-> leafblower
# =============================================================
cor_w <- cor(w_autumn, w_lb)
cat(sprintf("\n=== Weight agreement ===\nautumn ↔ leafblower Pearson r = %.4f\n", cor_w))

# =============================================================
# 4. Save reference for future use
# =============================================================
ref <- list(
  # Input metadata
  n             = length(w_autumn),
  K             = length(tgt),
  total_cats    = sum(margin_counts),
  margin_counts = margin_counts,
  margin_order  = names(tgt),
  data_path     = DATA_PATH,
  targets_path  = TARGETS_PATH,
  # Original Rmd parameters
  max_iter      = MAX_ITER,
  max_weight    = MAX_WEIGHT,
  min_weight    = 0,
  autumn_conv   = AUTUMN_CONV,
  # autumn reference (matches Rmd line 468-476 exactly)
  autumn_wall_s = autumn_wall_s,
  autumn_weights = w_autumn,
  autumn_deff   = autumn_deff,
  autumn_ess    = autumn_ess,
  autumn_max_w  = max(w_autumn),
  # leafblower comparison
  lb_wall_s     = lb_wall_s,
  lb_weights    = w_lb,
  lb_status     = info_lb$status,
  lb_iterations = info_lb$iterations,
  lb_max_error  = info_lb$max_error,
  lb_deff       = lb_deff,
  lb_ess        = lb_ess,
  # Agreement
  cor_autumn_lb = cor_w,
  # Provenance
  git_sha       = system("git rev-parse HEAD", intern=TRUE),
  R_version     = R.version.string,
  captured_at   = as.character(Sys.time())
)

out_rds <- file.path(OUT_DIR, "stepstone_reference.rds")
saveRDS(ref, out_rds)
cat(sprintf("\nReference saved: %s (%.1f MB)\n",
    out_rds, file.info(out_rds)$size / 1024^2))

# Compact summary (without full weight vectors) for quick inspection
summary_rds <- file.path(OUT_DIR, "stepstone_reference_summary.rds")
ref_summary <- ref
ref_summary$autumn_weights <- NULL
ref_summary$lb_weights     <- NULL
# Keep distribution summaries for quick lookup
ref_summary$autumn_w_summary <- summary(w_autumn)
ref_summary$lb_w_summary     <- summary(w_lb)
saveRDS(ref_summary, summary_rds)
cat(sprintf("Summary saved: %s\n", summary_rds))

cat("\nDone.\n")
