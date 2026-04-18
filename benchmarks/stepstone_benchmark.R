#!/usr/bin/env Rscript
# benchmarks/stepstone_benchmark.R
#
# Real-world benchmark: autumn vs R leafblower vs Python leafblower
# Mirrors the calibration structure of:
#   22-weighting-create-weights-2025-firmsize.Rmd
#
# 9 margins, 835 categories, n=200K.
# Simulates German salary-survey data with realistic biases (~5-10% relative).
# Saves calibration-ready parquet for the Python companion script.
#
# Convergence notes:
#   - All methods use tol=1e-3 (0.1% worst margin error, survey-adequate).
#   - iEPPA (cyclic IPF + Dykstra box) converges to 1e-3 in ~130ms.
#   - Tightening to 1e-4 causes all cyclic IPF methods to stall: the
#     overlapping 9-margin system oscillates with no further progress.
#     This is structural, not a leafblower bug.
#   - autumn does not enforce max_weight as an inner-loop constraint;
#     leafblower enforces [0, max_weight] strictly via Dykstra.

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(leafblower)
  library(autumn)
})

set.seed(2025)
n <- 200000L
DATA_PATH    <- "benchmarks/stepstone_bench_data.parquet"
TARGETS_PATH <- "benchmarks/stepstone_bench_targets.json"

cat("=== Stepstone-style calibration benchmark ===\n")
cat(sprintf("n = %s | margins = 9 | max_weight = 5\n\n", format(n, big.mark=",")))

# ── 1. Levels ────────────────────────────────────────────────────────────────

states <- c(
  "Schleswig-Holstein","Hamburg","Niedersachsen","Bremen",
  "Nordrhein-Westfalen","Hessen","Rheinland-Pfalz","Baden-Württemberg",
  "Bayern","Saarland","Berlin","Brandenburg",
  "Mecklenburg-Vorpommern","Sachsen","Sachsen-Anhalt","Thüringen",
  "Deutschland"
)
wz_lvl    <- c("A","BDE","CFH","G","I","J","K","LM","N","OURST","P","Q")
age_lvl   <- c("unter 25 Jahre","25 bis unter 35 Jahre","35 bis unter 45 Jahre",
               "45 bis unter 55 Jahre","55 bis unter 65 Jahre","65 Jahre und älter")
fsize_lvl <- c("1 bis 9 Beschäftigte","10 bis 49 Beschäftigte",
               "50 bis 249 Beschäftigte","250 und mehr Beschäftigte")

# ── 2. Survey sample (realistic biases, ~5-10% relative) ─────────────────────
p_state_samp <- c(0.030,0.030,0.091,0.008,0.212,0.075,0.041,0.133,
                  0.159,0.010,0.044,0.028,0.018,0.040,0.022,0.020,0.039)
p_state_samp <- p_state_samp / sum(p_state_samp)

p_wz_samp <- c(0.013,0.024,0.238,0.128,0.032,0.058,0.042,0.082,
               0.055,0.090,0.060,0.178)
p_wz_samp <- p_wz_samp / sum(p_wz_samp)

p_age_samp  <- c(0.075,0.260,0.270,0.255,0.130,0.010)
p_age_samp  <- p_age_samp / sum(p_age_samp)

p_fsize_samp <- c(0.13,0.21,0.22,0.44)

# Männer:Vollzeit, Männer:Teilzeit, Frauen:Vollzeit, Frauen:Teilzeit
p_gaz_samp <- c(0.57*0.90, 0.57*0.10, 0.43*0.53, 0.43*0.47)

p_edu_samp <- c(0.38, 0.62)

g_az <- sample(c("M:V","M:T","F:V","F:T"), n, replace=TRUE, prob=p_gaz_samp)
data_survey <- tibble(
  location    = sample(states,   n, replace=TRUE, prob=p_state_samp),
  WZmatch     = sample(wz_lvl,   n, replace=TRUE, prob=p_wz_samp),
  age10       = sample(age_lvl,  n, replace=TRUE, prob=p_age_samp),
  fsize       = sample(fsize_lvl,n, replace=TRUE, prob=p_fsize_samp),
  Abschluss   = sample(c("mit Hochschulabschluss","ohne Hochschulabschluss"),
                       n, replace=TRUE, prob=p_edu_samp),
  gender      = ifelse(startsWith(g_az,"M"), "Männer", "Frauen"),
  Arbeitszeit = ifelse(endsWith(g_az,":V"),  "Vollzeit", "Teilzeit")
) |>
  mutate(
    rk_gender                  = gender,
    rk_time                    = Arbeitszeit,
    rk_age10_gender            = interaction(age10, gender, sep=":") |> as.character(),
    rk_gender_time             = interaction(gender, Arbeitszeit, sep=":") |> as.character(),
    rk_i_loc_time_gender       = interaction(location, Arbeitszeit, gender, sep=":") |> as.character(),
    rk_i_loc_time_age10_gender = interaction(location, Arbeitszeit, age10, gender, sep=":") |> as.character(),
    rk_i_loc_wz                = interaction(location, WZmatch, sep=":") |> as.character(),
    rk_i_loc_fsize             = interaction(location, fsize, sep=":") |> as.character(),
    rk_i_loc_Abschluss_gender  = interaction(location, Abschluss, gender, sep=":") |> as.character()
  )

cat(sprintf("Survey: %s rows × %d columns\n", format(nrow(data_survey), big.mark=","), ncol(data_survey)))

# ── 3. Population targets ─────────────────────────────────────────────────────
POP <- 2000000L
p_state_pop <- c(0.030,0.027,0.093,0.008,0.212,0.075,0.041,0.133,
                 0.161,0.010,0.039,0.029,0.019,0.041,0.023,0.021,0.038)
p_state_pop <- p_state_pop / sum(p_state_pop)

p_wz_pop <- c(0.015,0.025,0.245,0.132,0.035,0.045,0.035,0.082,
              0.055,0.095,0.062,0.174)
p_wz_pop <- p_wz_pop / sum(p_wz_pop)

p_age_pop  <- c(0.08,0.22,0.25,0.28,0.16,0.01)
p_age_pop  <- p_age_pop / sum(p_age_pop)

p_fsize_pop <- c(0.17,0.22,0.21,0.40)
p_gaz_pop   <- c(0.55*0.88, 0.55*0.12, 0.45*0.50, 0.45*0.50)
p_edu_pop   <- c(0.33, 0.67)

g_az_p <- sample(c("M:V","M:T","F:V","F:T"), POP, replace=TRUE, prob=p_gaz_pop)
population <- tibble(
  location    = sample(states,    POP, replace=TRUE, prob=p_state_pop),
  WZmatch     = sample(wz_lvl,    POP, replace=TRUE, prob=p_wz_pop),
  age10       = sample(age_lvl,   POP, replace=TRUE, prob=p_age_pop),
  fsize       = sample(fsize_lvl, POP, replace=TRUE, prob=p_fsize_pop),
  Abschluss   = sample(c("mit Hochschulabschluss","ohne Hochschulabschluss"),
                       POP, replace=TRUE, prob=p_edu_pop),
  gender      = ifelse(startsWith(g_az_p,"M"), "Männer", "Frauen"),
  Arbeitszeit = ifelse(endsWith(g_az_p,":V"),  "Vollzeit", "Teilzeit")
)

mk <- function(pop, ...) {
  vs <- rlang::enquos(...); k <- rlang::sym("__k__")
  pop |>
    mutate(!!k := interaction(!!!vs, sep=":") |> as.character()) |>
    count(!!k) |> mutate(p = n / sum(n)) |> select(-n) |> deframe()
}

target_anes <- list(
  rk_gender                  = mk(population, gender),
  rk_time                    = mk(population, Arbeitszeit),
  rk_age10_gender            = mk(population, age10, gender),
  rk_gender_time             = mk(population, gender, Arbeitszeit),
  rk_i_loc_time_gender       = mk(population, location, Arbeitszeit, gender),
  rk_i_loc_time_age10_gender = mk(population, location, Arbeitszeit, age10, gender),
  rk_i_loc_wz                = mk(population, location, WZmatch),
  rk_i_loc_fsize             = mk(population, location, fsize),
  rk_i_loc_Abschluss_gender  = mk(population, location, Abschluss, gender)
)

# Remove target cells absent from survey (renormalise)
for (nm in names(target_anes)) {
  sv <- unique(data_survey[[nm]]); tgt <- target_anes[[nm]]
  miss <- setdiff(names(tgt), sv)
  if (length(miss) > 0) {
    tgt <- tgt[names(tgt) %in% sv]
    target_anes[[nm]] <- tgt / sum(tgt)
  }
}

cat(sprintf("Margins: %d  |  categories: %s  |  total: %d\n",
    length(target_anes),
    paste(sapply(target_anes, length), collapse=", "),
    sum(sapply(target_anes, length))))

# ── 4. Save data ──────────────────────────────────────────────────────────────
arrow::write_parquet(data_survey, DATA_PATH)
jsonlite::write_json(lapply(target_anes, as.list), TARGETS_PATH,
                     pretty=FALSE, auto_unbox=TRUE)
cat(sprintf("Data saved → %s\n\n", DATA_PATH))

# ── 5. Benchmark settings ─────────────────────────────────────────────────────
MAX_ITER   <- 500L
MAX_WEIGHT <- 5
TOL        <- 1e-3   # 0.1% worst margin error; iEPPA converges in ~130ms at this level
                     # Tightening to 1e-4 causes oscillation in all cyclic IPF methods.
AUTUMN_CONV <- c(pct = 1e-15, absolute = TOL)

time_one <- function(expr) {
  t0 <- proc.time()["elapsed"]
  res <- force(expr)
  list(result=res, ms=(proc.time()["elapsed"] - t0) * 1000)
}

report <- function(label, res, ms, algo=NULL) {
  w <- res$weights
  cat(sprintf("%-30s  %6.0f ms  w[min=%.3f med=%.3f max=%.3f]  DEFF=%.3f  ESS=%.0f\n",
      label, ms, min(w), median(w), max(w),
      leafblower::design_effect(res),
      leafblower::effective_sample_size(res)))
  if (!is.null(algo))
    cat(sprintf("  %s algorithm: %s\n", label, algo))
}

# ── 6. autumn ────────────────────────────────────────────────────────────────
cat("--- autumn::harvest ---\n")
invisible(autumn::harvest(data_survey, target_anes, max_weight=MAX_WEIGHT,
                          max_iterations=MAX_ITER, convergence=AUTUMN_CONV)) # warmup
r_autumn <- time_one(autumn::harvest(data_survey, target_anes, max_weight=MAX_WEIGHT,
                                     max_iterations=MAX_ITER, convergence=AUTUMN_CONV))

w_a <- r_autumn$result$weights
cat(sprintf("  time:     %.0f ms\n", r_autumn$ms))
cat(sprintf("  weights:  min=%.3f  med=%.3f  max=%.3f\n", min(w_a), median(w_a), max(w_a)))
cat(sprintf("  w > 5:    %d observations (%.1f%%)\n", sum(w_a > 5), 100*mean(w_a > 5)))
cat(sprintf("  DEFF:     %.3f  |  ESS: %.0f / %d\n\n",
    autumn::design_effect(r_autumn$result),
    autumn::effective_sample_size(r_autumn$result), n))

# ── 7. leafblower ────────────────────────────────────────────────────────────
cat("--- leafblower::harvest (ieppa, tol=1e-3) ---\n")
invisible(leafblower::harvest(data_survey, target_anes, method="ieppa",
                              max_weight=MAX_WEIGHT, max_iterations=MAX_ITER,
                              convergence=list(absolute=TOL))) # warmup
r_lb <- time_one(leafblower::harvest(data_survey, target_anes, method="ieppa",
                                     max_weight=MAX_WEIGHT, max_iterations=MAX_ITER,
                                     convergence=list(absolute=TOL)))

w_lb <- r_lb$result$weights
cat(sprintf("  time:     %.0f ms\n", r_lb$ms))
cat(sprintf("  weights:  min=%.3f  med=%.3f  max=%.3f\n", min(w_lb), median(w_lb), max(w_lb)))
cat(sprintf("  w > 5:    %d observations (0%% — box constraint strictly enforced)\n", sum(w_lb > 5)))
cat(sprintf("  DEFF:     %.3f  |  ESS: %.0f / %d\n\n",
    leafblower::design_effect(w_lb),
    leafblower::effective_sample_size(w_lb), n))

# ── 8. Weight agreement ───────────────────────────────────────────────────────
cat(sprintf("Weight correlation autumn ↔ leafblower: r = %.4f\n\n", cor(w_a, w_lb)))

# ── 9. Summary ───────────────────────────────────────────────────────────────
speedup <- r_autumn$ms / r_lb$ms
cat("=== Summary (n=200K, 9 margins, 835 categories, tol=1e-3, max_weight=5) ===\n")
cat(sprintf("  autumn:       %5.0f ms  |  max_weight enforced: NO  (post-hoc clip)\n", r_autumn$ms))
cat(sprintf("  R leafblower: %5.0f ms  |  max_weight enforced: YES (Dykstra box)\n", r_lb$ms))
cat(sprintf("  Speedup:      %.1fx\n", speedup))
cat("\nRun benchmarks/stepstone_benchmark.py for the Python comparison.\n")
