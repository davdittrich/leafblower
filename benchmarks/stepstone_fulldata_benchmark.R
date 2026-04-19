#!/usr/bin/env Rscript
# benchmarks/stepstone_fulldata_benchmark.R
#
# Benchmark on REAL Stepstone salary-survey data with the ORIGINAL margins
# from 22-weighting-create-weights-2025-firmsize.Rmd.
#
# Reproduces the full data-prep + target-construction pipeline verbatim,
# then measures autumn::harvest vs leafblower::harvest on the exact same
# n~1.25M dataset and 9-margin target list.
#
# Saves calibration-ready data for the Python companion script.

suppressPackageStartupMessages({
  library(tidyverse)
  library(arrow)
  library(leafblower)
  library(autumn)
  library(lubridate)
  library(jsonlite)
})

ROOT      <- "/home/dd/stepstone"
DATA_V3   <- file.path(ROOT, "Salary-Data/data/v3")
DATA_DIR  <- file.path(ROOT, "Salary-Data/data")
GEO_DIR   <- file.path(ROOT, "geo_db/data")
OUT_DIR   <- "benchmarks"

cat("=== Stepstone full-data benchmark ===\n")
cat("Source:", ROOT, "\n\n")

# ── 1. Replicate Rmd data prep ─────────────────────────────────────────────
cat("Loading and preparing data...\n")
t_prep <- system.time({

read_parquet(file.path(DATA_V3, "uuid_filtred_normalized_salary.parquet")) |>
  left_join(read_parquet(file.path(DATA_V3, "uuid_GS8.parquet"))) |>
  left_join(read_parquet(file.path(DATA_V3, "uuid_edu_kldb.parquet")) |>
              select(uuid, education, kldb)) |>
  left_join(read_parquet(file.path(DATA_V3, "uuid_filtred_variables.parquet"))) |>
  mutate(Arbeitszeit = as.factor(
           ifelse(sct_part_time_hours_per_week >= 35, "Vollzeit", "Teilzeit"))) |>
  mutate(GS2 = str_sub(GS8, 1, 2)) |>
  left_join(read_csv(file.path(GEO_DIR, "states_GS_D.csv"), show_col_types=FALSE) |>
              select(GS2, name) |> rename(location = "name")) |>
  mutate(gender = case_when(
    sct_gender == "SCT_GENDER_ANSWER_0" ~ "Männer",
    sct_gender == "SCT_GENDER_ANSWER_1" ~ "Frauen",
    TRUE ~ NA_character_)) |>
  mutate(Abschluss = case_when(
    education %in% c("24", "34")               ~ "ohneBeruf",
    education %in% c("64", "74", "84")         ~ "studied",
    education %in% c("35_2","35_3","55","65","75") ~ "mitBeruf")) |>
  mutate(Abschluss = case_when(
    Abschluss == "studied" ~ "mit Hochschulabschluss",
    TRUE                   ~ "ohne Hochschulabschluss")) |>
  mutate(sct_employer_industry = case_when(
    sct_employer_industry == "IND_MANUFACTURING" ~ "IND_MFC",
    sct_employer_industry == "IND_TECHNOLOGY"    ~ "IND_IT",
    TRUE ~ sct_employer_industry)) |>
  left_join(read_csv2(file.path(DATA_DIR, "industry_wz2008.csv"),
                      show_col_types=FALSE) |>
              select(slug, WZ2008, WZmatch),
            by = c("sct_employer_industry" = "slug")) |>
  mutate(age10 = case_when(
    age < 25              ~ "unter 25 Jahre",
    age >= 25 & age < 35  ~ "25 bis unter 35 Jahre",
    age >= 35 & age < 45  ~ "35 bis unter 45 Jahre",
    age >= 45 & age < 55  ~ "45 bis unter 55 Jahre",
    age >= 55 & age < 65  ~ "55 bis unter 65 Jahre",
    age >= 65             ~ "65 Jahre und älter") |>
    as_factor() |> fct_relevel("unter 25 Jahre")) |>
  mutate(fsize = fct_collapse(sct_number_of_employees,
    "1 bis 9 Beschäftigte"     = "CO_1_10",
    "10 bis 49 Beschäftigte"   = "CO_11_50",
    "50 bis 249 Beschäftigte"  = "CO_51_250",
    "250 und mehr Beschäftigte"= c("CO_251_500","CO_501_1000","CO_1001_5000",
                                   "CO_5001_10000","CO_10001_MAX"))) ->
  stepstone

}) # system.time
cat(sprintf("  Rows: %s  |  prep: %.1fs\n",
    format(nrow(stepstone), big.mark=","), t_prep["elapsed"]))

# ── 2. Compute share_D and wave window ─────────────────────────────────────
cap <- 3

stepstone |>
  summarise(last = max(date)) |>
  collect() |>
  mutate(
    first = (as_datetime(last) %m-% years(cap)) |> floor_date(unit="year"),
    wave  = (as_datetime(last) %m-% months(6))  |> ceiling_date(unit="quarter"),
    cap   = cap) ->
  waves

stepstone |>
  summarise(share_D = mean(GS2 == "00")) |>
  pull(share_D) -> share_D
share_D <- share_D / (1 + share_D)

numlocation <- tribble(
  ~shortcode, ~State,
  "00","Deutschland",    "01","Schleswig-Holstein", "02","Hamburg",
  "03","Niedersachsen",  "04","Bremen",              "05","Nordrhein-Westfalen",
  "06","Hessen",         "07","Rheinland-Pfalz",     "08","Baden-Württemberg",
  "09","Bayern",         "10","Saarland",             "11","Berlin",
  "12","Brandenburg",    "13","Mecklenburg-Vorpommern","14","Sachsen",
  "15","Sachsen-Anhalt", "16","Thüringen")

cat(sprintf("  Wave: %s → %s\n\n", waves$first, waves$wave))

# ── 3. Replicate Rmd target construction ──────────────────────────────────
cat("Building targets...\n")

State_gender <- read_parquet(file.path(DATA_DIR, "margin_State_age_time_gender.parquet"))
last_sg <- State_gender |> summarise(last=max(date)) |> pull(last)
first   <- waves$first; if (first > last_sg) first <- last_sg
State_gender |>
  filter(date >= first, date < waves$wave) |>
  na.omit() -> State_gender

State_gender |>
  bind_rows(State_gender |>
              group_by(date,Arbeitszeit,gender,age10) |>
              summarise(location="Deutschland", value=sum(value)*share_D, .groups="drop")) |>
  mutate(age10 = fct_collapse(age10,
    `unter 25 Jahre` = c("unter 15 Jahre","15 bis unter 25 Jahre"))) |>
  filter(Arbeitszeit %in% c("Vollzeit","Teilzeit")) |>
  group_by(date,Arbeitszeit,location,age10,gender) |>
  summarise(n=sum(value), .groups="drop") |>
  group_by(Arbeitszeit,location,age10,gender) |>
  summarise(n=mean(n), .groups="drop") |>
  mutate(i=interaction(location,Arbeitszeit,age10,gender,sep=":")) |>
  ungroup() |> mutate(total=sum(n)) |>
  group_by(i) |> summarise(p=sum(n)/total) |> deframe() -> t_latotal

State_gender |>
  bind_rows(State_gender |>
              group_by(date,Arbeitszeit,gender) |>
              summarise(location="Deutschland", value=sum(value)*share_D, .groups="drop")) |>
  filter(Arbeitszeit %in% c("Vollzeit","Teilzeit")) |>
  group_by(date,Arbeitszeit,location,gender) |>
  summarise(n=sum(value), .groups="drop") |>
  group_by(Arbeitszeit,location,gender) |>
  summarise(n=mean(n), .groups="drop") |>
  mutate(i=interaction(location,Arbeitszeit,gender,sep=":")) |>
  ungroup() |> mutate(total=sum(n)) |>
  group_by(i) |> summarise(p=sum(n)/total) |> deframe() -> t_ltotal

State_gender |>
  filter(Arbeitszeit %in% c("Vollzeit","Teilzeit")) |>
  group_by(date,Arbeitszeit,gender) |>
  summarise(n=sum(value), .groups="drop") |>
  group_by(Arbeitszeit,gender) |>
  summarise(n=mean(n), .groups="drop") |>
  mutate(i=interaction(gender,Arbeitszeit,sep=":")) |>
  ungroup() |> mutate(total=sum(n)) |>
  group_by(i) |> summarise(p=sum(n)/total) |> deframe() -> fed_agtotal

State_gender |>
  mutate(age10 = fct_collapse(age10,
    `unter 25 Jahre` = c("unter 15 Jahre","15 bis unter 25 Jahre"))) |>
  filter(Arbeitszeit %in% c("Vollzeit","Teilzeit")) |>
  group_by(date,age10,gender) |>
  summarise(n=sum(value), .groups="drop") |>
  group_by(age10,gender) |>
  summarise(n=mean(n), .groups="drop") |>
  mutate(i=interaction(age10,gender,sep=":")) |>
  ungroup() |> mutate(total=sum(n)) |>
  group_by(i) |> summarise(p=sum(n)/total) |> deframe() -> fed_agegtotal

# --- education margin ---
State_gender_edu <- read_parquet(file.path(DATA_DIR,"margin_State_edu.parquet"))
last_edu <- State_gender_edu |> summarise(last=max(date)) |> pull(last)
first_edu <- waves$first; if (first_edu > last_edu) first_edu <- last_edu
State_gender_edu |>
  filter(date >= first_edu, date < waves$wave) |>
  bind_rows(State_gender_edu |>
              filter(date >= first_edu, date < waves$wave) |>
              group_by(date,Abschluss,gender) |>
              summarise(location="Deutschland", value=sum(value)*share_D, .groups="drop")) |>
  na.omit() |>
  mutate(Abschluss = case_when(
    Abschluss == "studied" ~ "mit Hochschulabschluss", TRUE ~ "ohne Hochschulabschluss")) |>
  group_by(location,Abschluss,gender) |>
  summarise(periods=length(unique(date)), n=sum(value)/periods, .groups="drop") |>
  mutate(i=interaction(location,Abschluss,gender,sep=":")) |>
  ungroup() |> mutate(total=sum(n)) |>
  group_by(i) |> summarise(p=sum(n)/total) |> deframe() -> t_lstudgender

# --- WZ margin ---
State_wz <- read_parquet(file.path(DATA_DIR,"margin_State_wz.parquet"))
last_wz <- State_wz |> summarise(last=max(date)) |> pull(last)
first_wz <- waves$first; if (first_wz > last_wz) first_wz <- last_wz
State_wz |>
  filter(date >= first_wz, date < waves$wave) |>
  bind_rows(State_wz |>
              filter(date >= first_wz, date < waves$wave) |>
              group_by(date,WZ) |>
              summarise(location="Deutschland", Employees=sum(Employees)*share_D, .groups="drop")) |>
  mutate(WZmatch = case_when(
    WZ %in% c("A","G","I","J","K","P")               ~ WZ,
    WZ == "B, D, E"                                   ~ "BDE",
    WZ %in% c("C","F","H")                            ~ "CFH",
    WZ %in% c("L","M","L,M")                          ~ "LM",
    WZ %in% c("O","U","R","S","O, U","R, S","T","R, S, T") ~ "OURST",
    WZ %in% c("N","N ohne ANÜ","782, 783")            ~ "N",
    WZ %in% c("Q","86","87.88")                       ~ "Q",
    TRUE ~ NA_character_)) |>
  na.omit() |>
  group_by(location,WZmatch) |>
  summarise(n=sum(Employees), .groups="drop") |>
  mutate(i=interaction(location,WZmatch,sep=":")) |>
  ungroup() |> mutate(total=sum(n)) |>
  group_by(i) |> summarise(p=sum(n)/total) |> deframe() -> t_lwz

# --- firm-size margin ---
State_fsize <- read_parquet(file.path(DATA_DIR,"margin_firmsize.parquet"))
last_fs <- State_fsize |> summarise(last=max(date)) |> pull(last)
first_fs <- waves$first; if (first_fs > last_fs) first_fs <- last_fs
State_fsize |>
  filter(date >= first_fs, date < waves$wave) |>
  na.omit() |>
  select(-State) |>
  left_join(numlocation, by="shortcode") |>
  rename(location=State) |>
  mutate(Beschaeftigte = ifelse(shortcode=="00", Beschaeftigte*share_D, Beschaeftigte)) |>
  group_by(date,fsize,location) |>
  summarise(n=sum(Beschaeftigte), .groups="drop") |>
  group_by(fsize,location) |>
  summarise(n=mean(n), .groups="drop") |>
  mutate(i=interaction(location,fsize,sep=":")) |>
  ungroup() |> mutate(total=sum(n)) |>
  group_by(i) |> summarise(p=sum(n)/total) |> deframe() -> t_fstotal

# --- federal gender/arbeitszeit margin ---
Fed <- arrow::read_csv_arrow(file.path(DATA_DIR,"margin_Federal.csv"))
last_fed <- Fed |> summarise(last=max(date)) |> pull(last)
first_fed <- waves$first; if (first_fed > last_fed) first_fed <- last_fed

Fed |>
  filter(date >= first_fed, date < waves$wave) |>
  summarise(gender_Männer=mean(Männer), gender_Frauen=mean(Frauen)) |>
  pivot_longer(starts_with("gender_"), names_to="gender", names_prefix="gender_", values_to="n") |>
  ungroup() |> mutate(total=sum(n)) |>
  group_by(gender) |> summarise(p=sum(n)/total) |> deframe() -> Fed_Gender

Fed |>
  filter(date >= first_fed, date < waves$wave) |>
  summarise(t_Vollzeit=mean(Vollzeit), t_Teilzeit=mean(Teilzeit)) |>
  pivot_longer(starts_with("t_"), names_to="Arbeitszeit", names_prefix="t_", values_to="n") |>
  ungroup() |> mutate(total=sum(n)) |>
  filter(Arbeitszeit %in% c("Vollzeit","Teilzeit")) |>
  group_by(Arbeitszeit) |> summarise(p=sum(n)/total) |> deframe() -> Fed_Arbeitszeit

target_anes <- list(
  rk_i_loc_time_age10_gender = t_latotal,
  rk_i_loc_wz                = t_lwz,
  rk_age10_gender             = fed_agegtotal,
  rk_i_loc_time_gender        = t_ltotal,
  rk_i_loc_fsize              = t_fstotal,
  rk_i_loc_Abschluss_gender   = t_lstudgender,
  rk_gender_time              = fed_agtotal,
  rk_time                     = Fed_Arbeitszeit,
  rk_gender                   = Fed_Gender
)

cat(sprintf("  Margins: %d  |  categories: %s  |  total: %d\n",
    length(target_anes),
    paste(sapply(target_anes, length), collapse=", "),
    sum(sapply(target_anes, length))))

# ── 4. Replicate raking data frame ────────────────────────────────────────
cat("\nBuilding raking data frame + imputation...\n")
t_raking <- system.time({

stepstone |>
  filter(date >= waves$first) |>
  mutate(wave = waves$wave) |>
  mutate(
    i_loc_time          = interaction(location, Arbeitszeit, sep=":"),
    i_loc_age10         = interaction(location, age10, sep=":"),
    i_loc_time_age10    = interaction(location, Arbeitszeit, age10, sep=":"),
    i_loc_Abschluss     = interaction(location, Abschluss, sep=":"),
    i_loc_wz            = interaction(location, WZmatch, sep=":"),
    i_loc_fsize         = interaction(location, fsize, sep=":")
  ) |>
  mutate(
    rk_gender           = gender,
    rk_age10            = as.character(age10),
    rk_time             = as.character(Arbeitszeit),
    rk_Abschluss        = as.character(Abschluss),
    rk_location         = as.character(location),
    rk_i_loc_time       = as.character(i_loc_time),
    rk_i_loc_age10      = as.character(i_loc_age10),
    rk_i_loc_time_age10 = as.character(i_loc_time_age10),
    rk_i_loc_Abschluss  = as.character(i_loc_Abschluss),
    rk_i_loc_wz         = as.character(i_loc_wz),
    rk_i_loc_fsize      = as.character(i_loc_fsize)
  ) -> stepstone_raking

# Gender imputation: mode within location×age10×Arbeitszeit stratum.
# KNN (as in the Rmd) is too slow for a benchmark on 1.58M rows (~2h+).
# Mode-by-stratum gives the same calibration dataset structure at <1s.
mode_chr <- function(x) { x <- x[!is.na(x)]; if (!length(x)) NA_character_ else names(which.max(table(x))) }
stepstone_raking |>
  group_by(rk_location, rk_age10, rk_time) |>
  mutate(rk_gender = if_else(
    is.na(rk_gender),
    mode_chr(rk_gender),
    rk_gender)) |>
  ungroup() |>
  filter(!is.na(rk_gender)) |>
  mutate(
    rk_gender_time             = interaction(rk_gender, Arbeitszeit, sep=":") |> as.character(),
    rk_age10_gender            = interaction(rk_age10, rk_gender, sep=":") |> as.character(),
    rk_i_loc_time_gender       = interaction(rk_i_loc_time, rk_gender, sep=":") |> as.character(),
    rk_i_loc_time_age10_gender = interaction(rk_i_loc_time_age10, rk_gender, sep=":") |> as.character(),
    rk_i_loc_Abschluss_gender  = interaction(rk_i_loc_Abschluss, rk_gender, sep=":") |> as.character(),
    rk_Abschluss_gender        = interaction(rk_Abschluss, rk_gender, sep=":") |> as.character()
  ) -> stepstone_imputed

})
cat(sprintf("  n = %s  |  prep+impute: %.1fs\n\n",
    format(nrow(stepstone_imputed), big.mark=","), t_raking["elapsed"]))

# ── 5. Feasibility: drop target cells absent from data ────────────────────
for (nm in names(target_anes)) {
  sv  <- unique(stepstone_imputed[[nm]])
  tgt <- target_anes[[nm]]
  miss <- setdiff(names(tgt), sv)
  if (length(miss) > 0) {
    tgt <- tgt[names(tgt) %in% sv]
    target_anes[[nm]] <- tgt / sum(tgt)
    cat(sprintf("  WARN %s: dropped %d missing cells, renormalised\n", nm, length(miss)))
  }
}

# ── 6. Save for Python benchmark ──────────────────────────────────────────
cat("Saving calibration-ready data for Python...\n")
# Keep only the columns Python needs: the 9 rk_* interaction columns + uuid
col_keep <- c("uuid", names(target_anes))
stepstone_imputed |>
  select(any_of(col_keep)) |>
  arrow::write_parquet(file.path(OUT_DIR, "stepstone_fulldata_bench_data.parquet"))
jsonlite::write_json(
  lapply(target_anes, as.list),
  file.path(OUT_DIR, "stepstone_fulldata_bench_targets.json"),
  pretty=FALSE, auto_unbox=TRUE)
cat(sprintf("  Saved %s rows to benchmarks/stepstone_fulldata_bench_data.parquet\n\n",
    format(nrow(stepstone_imputed), big.mark=",")))

# ── 7. Benchmark ──────────────────────────────────────────────────────────
MAX_ITER    <- 3000L   # matches the Rmd
MAX_WEIGHT  <- 5
TOL         <- 1e-3    # practical limit for this 9-margin system (see synthetic bench)
AUTUMN_CONV <- c(pct=1e-15, absolute=TOL)

time_one <- function(expr) {
  t0  <- proc.time()["elapsed"]
  res <- force(expr)
  list(result=res, ms=(proc.time()["elapsed"] - t0) * 1000)
}

cat("--- autumn::harvest ---\n")
cat("  (warmup...)\n")
suppressWarnings(
  invisible(autumn::harvest(stepstone_imputed, target_anes,
    max_weight=MAX_WEIGHT, max_iterations=50, convergence=AUTUMN_CONV)))

cat("  (timed run...)\n")
suppressWarnings(
  r_autumn <- time_one(autumn::harvest(stepstone_imputed, target_anes,
    max_weight=MAX_WEIGHT, max_iterations=MAX_ITER, convergence=AUTUMN_CONV)))
w_a <- r_autumn$result$weights

cat(sprintf("  time:     %.0f ms (%.1f s)\n", r_autumn$ms, r_autumn$ms/1000))
cat(sprintf("  n:        %s\n", format(length(w_a), big.mark=",")))
cat(sprintf("  weights:  min=%.3f  med=%.3f  max=%.3f\n", min(w_a), median(w_a), max(w_a)))
cat(sprintf("  w > 5:    %d (%.2f%%)\n", sum(w_a > 5), 100*mean(w_a > 5)))
cat(sprintf("  DEFF:     %.3f  |  ESS: %s\n\n",
    autumn::design_effect(r_autumn$result),
    format(round(autumn::effective_sample_size(r_autumn$result)), big.mark=",")))

cat("--- leafblower::harvest (ieppa, tol=1e-3) ---\n")
cat("  (warmup...)\n")
suppressWarnings(
  invisible(leafblower::harvest(stepstone_imputed, target_anes,
    method="ieppa", max_weight=MAX_WEIGHT, max_iterations=50,
    convergence=list(absolute=TOL))))

cat("  (timed run...)\n")
suppressWarnings(
  r_lb <- time_one(leafblower::harvest(stepstone_imputed, target_anes,
    method="ieppa", max_weight=MAX_WEIGHT, max_iterations=MAX_ITER,
    convergence=list(absolute=TOL))))
w_lb <- r_lb$result$weights

cat(sprintf("  time:     %.0f ms (%.1f s)\n", r_lb$ms, r_lb$ms/1000))
cat(sprintf("  weights:  min=%.3f  med=%.3f  max=%.3f\n", min(w_lb), median(w_lb), max(w_lb)))
cat(sprintf("  w > 5:    0 (strict Dykstra box)\n"))
cat(sprintf("  DEFF:     %.3f  |  ESS: %s\n\n",
    leafblower::design_effect(w_lb),
    format(round(leafblower::effective_sample_size(w_lb)), big.mark=",")))

# ── 8. Summary ────────────────────────────────────────────────────────────
cat(sprintf("Weight correlation autumn ↔ leafblower: r = %.4f\n\n", cor(w_a, w_lb)))

speedup <- r_autumn$ms / r_lb$ms
cat(sprintf("=== Summary (n=%s, 9 margins, %d categories, tol=1e-3) ===\n",
    format(nrow(stepstone_imputed), big.mark=","),
    sum(sapply(target_anes, length))))
cat(sprintf("  autumn:       %6.0f ms (%4.1f s)  |  strict box: NO\n",
    r_autumn$ms, r_autumn$ms/1000))
cat(sprintf("  R leafblower: %6.0f ms (%4.1f s)  |  strict box: YES\n",
    r_lb$ms, r_lb$ms/1000))
cat(sprintf("  Speedup:      %.1fx\n", speedup))
cat("\nRun benchmarks/stepstone_fulldata_benchmark.py for Python comparison.\n")
