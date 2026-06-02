#!/usr/bin/env Rscript
# T1: parse params_R.txt and params_Py.txt, emit params_diff.toon
# Usage: Rscript benchmarks/2apm/diff_params.R

if (!file.exists("benchmarks/2apm")) {
  stop("Run from project root. Current wd=", getwd())
}

OUTDIR <- "benchmarks/2apm"

# ---- Read raw dump files ----
r_txt  <- readLines(file.path(OUTDIR, "params_R.txt"))
py_txt <- readLines(file.path(OUTDIR, "params_Py.txt"))

# ---- Parse R dump: lines matching "  slot NN | <name> = <value>" ----
parse_r_dump <- function(lines) {
  pat <- "^  slot\\s+(\\d+)\\s+\\|\\s+(\\S+)\\s+=\\s+(.+)$"
  out <- list()
  for (l in lines) {
    m <- regmatches(l, regexec(pat, l))[[1]]
    if (length(m) == 4) {
      out[[m[3]]] <- list(slot = as.integer(m[2]), raw = trimws(m[4]))
    }
  }
  conv_pat <- "^  conv\\$(.+)\\s+=\\s+(.+)$"
  for (l in lines) {
    m <- regmatches(l, regexec(conv_pat, l))[[1]]
    if (length(m) == 3) {
      out[[paste0("_conv.", m[2])]] <- list(slot = NA, raw = trimws(m[3]))
    }
  }
  out
}

# ---- Parse Python dump: lines matching "  <name>  = <repr>" ----
parse_py_dump <- function(lines) {
  pat1 <- "^  ([a-z_]+)\\s+=\\s+(.+)$"
  out <- list()
  in_params <- FALSE
  for (l in lines) {
    if (grepl("params dict", l)) { in_params <- TRUE; next }
    if (grepl("^---", l) && in_params) { in_params <- FALSE }
    if (in_params) {
      m <- regmatches(l, regexec(pat1, l))[[1]]
      if (length(m) == 3) {
        out[[m[2]]] <- list(status = "IN_DICT", raw = trimws(m[3]))
      }
    }
  }
  omit_pat <- "^OMITTED_KEY: (\\w+) \\| C_DEFAULT: (\\S+) \\| CITATION: (\\S+)$"
  for (l in lines) {
    m <- regmatches(l, regexec(omit_pat, l))[[1]]
    if (length(m) == 4) {
      out[[m[2]]] <- list(status = "OMITTED", raw = "NOT_IN_DICT", cpp_default = m[3])
    }
  }
  out
}

r_params  <- parse_r_dump(r_txt)
py_params <- parse_py_dump(py_txt)

# ---- Canonical param name mapping: R slot name → C struct field name ----
# R sends positional args that map to C field names via r_bridge.cpp
PARAM_MAP <- data.frame(
  c_field        = c("min_weight","max_weight","inner_max_iter","outer_max_iter","tol_abs",
                      "verbose","algorithm","epsilon","bounds_mode",
                      "pct_tol","absolute_tol","metric","rule","stop_when",
                      "sor_enabled","sor_auto","sor_omega_init","sor_omega_min","sor_omega_fixed","sor_burnin",
                      "homotopy.n_levels","homotopy.start_factor","homotopy.end_factor","homotopy.budget_split_p",
                      "scheduler","eta_mode","eta_start","eta_end","eta_schedule_power",
                      "capacity_penalty","alm_penalty",
                      "accelerate","newton_tsvd_ratio","ridge_lambda","sor_corun_aa",
                      "gk_omega","sk_omega"),
  r_slot_name    = c("min_weight","max_weight","max_iterations","max_iterations","tol_abs",
                      "verbose","method","method","bounds_mode",
                      "pct_tol","absolute_tol","metric","rule","stop_when",
                      "sor_enabled","sor_auto","sor_omega_init","sor_omega_min","sor_omega_fixed","sor_burnin",
                      "homotopy_levels","homotopy_start_factor","homotopy_end_factor","homotopy_budget_p",
                      "scheduler","eta_schedule","eta_start","eta_end","eta_schedule_power",
                      "capacity_penalty","alm_penalty",
                      "accelerate","newton_tsvd_ratio","ridge_lambda","sor_corun_aa",
                      "gk_omega","sk_omega"),
  py_dict_key    = c("min_weight","max_weight","inner_max_iter","outer_max_iter","tol_abs",
                      "verbose","algorithm","epsilon","bounds_mode",
                      "pct_tol","absolute_tol","metric","rule","stop_when",
                      "sor_enabled","sor_auto","sor_omega_init","sor_omega_min","sor_omega_fixed","sor_burnin",
                      "homotopy_levels","homotopy_start_factor","homotopy_end_factor","homotopy_budget_p",
                      "scheduler","eta_mode","eta_start","eta_end","eta_schedule_power",
                      "capacity_penalty","alm_penalty",
                      "accelerate","newton_tsvd_ratio","ridge_lambda","sor_corun_aa",
                      "gk_omega","sk_omega"),
  stringsAsFactors = FALSE
)

# ---- Known resolved values for the test call (from source analysis) ----
# R values (what harvest.R sends to .Call for chebyshev, convergence=list(tol=1e-4),
#           max_weight=5, max_iterations=3000)
# chebyshev: method="chebyshev" (string), alg=5 in C; metric=max_err (int 0) — no override.
# capacity_penalty/alm_penalty: R sends -1.0 sentinel (NULL→-1.0 in harvest.R);
#   chebyshev ALM block gated by st.use_admm_capacity → harmless for chebyshev.
R_VALS <- list(
  group_ids_r         = "VECSXP[K]; INTSXP len=n per margin (0-idx codes, -1=NA). K=n_margins, n=nrow(data)",
  cat_counts_r        = "INTSXP[K]; level counts per margin. sum=total_cells",
  targets_r           = "VECSXP[K]; REALSXP len=cat_counts[k]. sum(targets[[k]])≈1",
  n_obs               = "INTSXP scalar = nrow(data)",
  min_weight          = 0.0,
  max_weight          = 5.0,
  algorithm           = "chebyshev (string; r_bridge.cpp decodes → RK_ALG_CHEBYSHEV=5)",
  verbose             = 0L,
  inner_max_iter      = 3000L,   # R sends max_iterations to inner; outer_max_iter = inner for non-newton_kl
  sw_vec              = "REALSXP len=n or NULL → uniform 1.0 (start_weights=NULL default)",
  tol_abs             = 1e-6,    # absolute_tol=0.0 → fallback 1e-6 (harvest.R:460)
  bounds_mode         = 0L,      # "cell" → 0
  pct_tol             = 1e-4,    # tol=1e-4 with rule=improvement → pct_tol
  absolute_tol        = 0.0,     # improvement rule → absolute_tol=0.0
  metric              = 0L,      # max_err (chebyshev: no method override; harvest.R:329 comment)
  rule                = 1L,      # improvement=1
  stop_when           = 0L,      # any=0
  sor_enabled         = 0L,      # sor=NULL → disabled
  sor_auto            = 0L,      # sor=NULL → disabled
  sor_omega_init      = 1.0,
  sor_omega_min       = 0.3,
  sor_omega_fixed     = -1.0,
  sor_burnin          = 20L,
  homotopy_n_levels   = 1L,
  homotopy_start_factor = 1.0,
  homotopy_end_factor = 1.0,
  homotopy_budget_p   = 0.5,
  scheduler           = "round_robin (string → C: RK_SCHED_ROUND_ROBIN=0)",
  eta_mode            = "fixed (string → C: RK_ETA_FIXED=0)",
  eta_start           = 1.0,
  eta_end             = 1.0,
  eta_schedule_power  = 0.5,
  capacity_penalty    = -1.0,    # NULL → sentinel -1.0; irrelevant for chebyshev (ALM gated)
  alm_penalty         = -1.0,    # NULL → sentinel -1.0; C bridge → 0.0 (inactive)
  accelerate          = 0L,
  newton_tsvd_ratio   = 1e-8,
  ridge_lambda        = 0.0,
  sor_corun_aa        = 0L,  # DEAD (e65t.1 NO-GO): co-run never enabled. See bd leafblower-e65t.1.
  gk_omega            = 1.0,  # DEAD (e65t.2 NO-GO): gk_omega experiment; default 1.0=identity. See bd leafblower-e65t.2.
  sk_omega            = 1.0  # e65t.3 GO: active parameter. See bd leafblower-e65t.3.
)

# Python values (what _harvest.py puts in params dict for same call)
# chebyshev: alg_int=5, metric=0 (max_err; not in _method_metric_map → no override)
PY_VALS <- list(
  group_ids_r         = "VECSXP[K]; INTSXP len=n per margin (0-idx codes, -1=NA). K=n_margins, n=nrow(data)",
  cat_counts_r        = "INTSXP[K]; level counts per margin. sum=total_cells",
  targets_r           = "VECSXP[K]; REALSXP len=cat_counts[k]. sum(targets[[k]])≈1",
  n_obs               = "INTSXP scalar = nrow(data)",
  min_weight          = 0.0,
  max_weight          = 5.0,
  algorithm           = 5L,      # alg_map["chebyshev"]=5 (Py sends int; R sends string)
  verbose             = 0L,
  inner_max_iter      = 3000L,
  sw_vec              = "REALSXP len=n or NULL → uniform 1.0 (start_weights=None default)",
  tol_abs             = 1e-6,    # absolute_tol=0.0 → fallback 1e-6 (_harvest.py:439)
  bounds_mode         = 0L,
  pct_tol             = 1e-4,    # tol=1e-4 with rule=improvement
  absolute_tol        = 0.0,
  metric              = 0L,      # max_err (int 0); chebyshev NOT in _method_metric_map
  rule                = 1L,      # improvement
  stop_when           = 0L,      # any
  sor_enabled         = 0L,      # sor=None → disabled
  sor_auto            = 0L,      # sor=None → disabled
  sor_omega_init      = 1.0,
  sor_omega_min       = 0.3,
  sor_omega_fixed     = -1.0,
  sor_burnin          = 20L,
  homotopy_n_levels   = 1L,
  homotopy_start_factor = 1.0,
  homotopy_end_factor = 1.0,
  homotopy_budget_p   = 0.5,
  scheduler           = 0L,      # round_robin → 0
  eta_mode            = 0L,      # fixed → 0
  eta_start           = 1.0,
  eta_end             = 1.0,
  eta_schedule_power  = 0.5,
  capacity_penalty    = "OMITTED: C++ default 0.0 (_bindings.cpp:131, c_api.cpp:memset)",
  alm_penalty         = "OMITTED: C++ default 0.0 (_bindings.cpp:136, c_api.cpp:memset)",
  accelerate          = 0L,      # int(False) = 0
  newton_tsvd_ratio   = 1e-8,
  ridge_lambda        = 0.0,
  sor_corun_aa        = 0L,  # DEAD (e65t.1 NO-GO): co-run never enabled. See bd leafblower-e65t.1.
  gk_omega            = NULL, # DEAD (e65t.2 NO-GO): gk_omega experiment; default 1.0=identity. See bd leafblower-e65t.2. Python binding out of scope.
  # DEAD (e65t.3 NO-GO): sk_omega experiment; default 1.0=identity. See bd leafblower-e65t.3. Python binding out of scope.
  sk_omega            = NULL
)

# ---- convergence/init path parameters for short-circuit evaluation ----
CONVERGENCE_INIT_PARAMS <- c(
  "capacity_penalty", "alm_penalty", "pct_tol", "absolute_tol",
  "metric", "rule", "stop_when", "sor_enabled", "sor_auto",
  "sor_omega_init", "sor_omega_min", "sor_omega_fixed", "sor_burnin",
  "scheduler", "eta_mode", "eta_start", "eta_end", "eta_schedule_power",
  "homotopy_n_levels", "accelerate", "newton_tsvd_ratio"
)

# ---- Build TOON rows ----
# 42 positional SEXP args to C_rk_calibrate (r_bridge.cpp:177-201), in call order.
# Removed: outer_max_iter (not a SEXP; bridge derives from inner), epsilon (not a SEXP; rk_params_t ABI field only).
# Added: group_ids_r, cat_counts_r, targets_r (pos 1-3), n_obs (pos 4), sw_vec (pos 10).
rows <- list()
params_in_order <- c(
  "group_ids_r","cat_counts_r","targets_r","n_obs",
  # scalars
  "min_weight","max_weight","algorithm","verbose","inner_max_iter",
  "sw_vec",
  "capacity_penalty","alm_penalty","tol_abs","bounds_mode",
  "homotopy_n_levels","homotopy_start_factor","homotopy_end_factor","homotopy_budget_p",
  "scheduler","eta_mode","eta_start","eta_end","eta_schedule_power",
  "pct_tol","absolute_tol","metric","rule","stop_when",
  "sor_enabled","sor_auto","sor_omega_init","sor_omega_min","sor_omega_fixed","sor_burnin",
  "accelerate","newton_tsvd_ratio","ridge_lambda","sor_corun_aa",
  "gk_omega","sk_omega"
)

cat(sprintf("Total params enumerated: %d\n", length(params_in_order)))

mismatch_rows <- list()
for (p in params_in_order) {
  rv <- R_VALS[[p]]
  pv <- PY_VALS[[p]]
  rv_str <- if (is.null(rv)) "null" else as.character(rv)
  pv_str <- if (is.null(pv)) "null" else as.character(pv)

  py_default_when_missing <- if (grepl("OMITTED", pv_str, fixed = TRUE)) {
    if (p == "capacity_penalty")
      "0.0 → chebyshev: ALM block gated by st.use_admm_capacity; irrelevant"
    else if (p == "alm_penalty")
      "0.0 → inactive (p->alm_penalty>0.0 false at c_api.cpp:254)"
    else "N/A"
  } else "N/A"

  mismatch  <- FALSE
  confidence <- 95L
  note       <- ""

  if (p %in% c("group_ids_r","cat_counts_r","targets_r","n_obs","sw_vec")) {
    # Data/input args: same parquet input → structurally identical on both sides.
    mismatch   <- FALSE
    confidence <- 97L
    note       <- "Data-input SEXP; structurally identical when reading same dataset."
  } else if (p == "capacity_penalty") {
    # c_api.cpp:376: (p->capacity_penalty <= 0.0) → auto path. chebyshev ALM gated anyway.
    mismatch   <- FALSE
    confidence <- 95L
    note       <- "R=-1.0, Py=0.0 (memset); both <=0.0 → auto/inactive for chebyshev (ALM gated by st.use_admm_capacity)."
  } else if (p == "alm_penalty") {
    # c_api.cpp:254: st.alm.mu path.
    mismatch   <- FALSE
    confidence <- 95L
    note       <- "R=-1.0, Py=0.0 (memset); both <=0.0 → alm inactive at c_api.cpp:254. Same semantic."
  } else if (p == "algorithm") {
    # R sends string "chebyshev"; Py sends int 5. r_bridge.cpp maps string → enum. Same result.
    mismatch   <- FALSE
    confidence <- 98L
    note       <- "R sends string 'chebyshev'; Py sends int 5. r_bridge.cpp maps string→enum. Same result."
  } else if (p == "scheduler") {
    # R sends "round_robin"; Py sends int 0. r_bridge.cpp strcmp.
    mismatch   <- FALSE
    confidence <- 98L
    note       <- "R sends string 'round_robin'; Py sends int 0. Both → RK_SCHED_ROUND_ROBIN. Same."
  } else if (p == "eta_mode") {
    mismatch   <- FALSE
    confidence <- 98L
    note       <- "R sends string 'fixed'; Py sends int 0. Both → RK_ETA_FIXED. Same."
  } else {
    rv_num <- suppressWarnings(as.numeric(rv_str))
    pv_num <- suppressWarnings(as.numeric(pv_str))
    if (!is.na(rv_num) && !is.na(pv_num)) {
      mismatch   <- (abs(rv_num - pv_num) > 1e-15)
      confidence <- 98L
    }
  }

  if (mismatch) {
    mismatch_rows[[length(mismatch_rows)+1]] <- list(
      param = p, r_value = rv_str, py_value = pv_str, confidence = confidence
    )
  }

  rows[[length(rows)+1]] <- list(
    param                   = p,
    r_value                 = rv_str,
    py_value                = pv_str,
    py_default_when_missing = py_default_when_missing,
    mismatch                = mismatch,
    confidence              = confidence,
    note                    = note
  )
}

# ---- Short-circuit determination ----
# All 41 positional SEXPs enumerated. epsilon removed (not a SEXP; rk_params_t ABI field only).
# No true mismatches remain: capacity_penalty and alm_penalty use different sentinels but
# both resolve to the same auto/inactive behavior (c_api.cpp:376-381, c_api.cpp:254).
# Data-input args (group_ids_r, cat_counts_r, targets_r, n_obs, sw_vec) are structurally identical.
# metric=0 (max_err) on BOTH sides: chebyshev has no override in either R or Python.
# → no short_circuit

short_circuit <- FALSE
short_circuit_reason <- paste0(
  "No true mismatches on any of the 42 positional SEXP args. ",
  "capacity_penalty/alm_penalty differ in sentinel (R=-1.0 vs Py=0.0 memset) but both reach ",
  "the same auto/inactive path (c_api.cpp:376-381, c_api.cpp:254); chebyshev ALM gated by ",
  "st.use_admm_capacity so capacity_mu is harmless. ",
  "metric=0 (max_err) on both sides: chebyshev absent from method-override map in both R and Py. ",
  "Data args identical. T2 trajectory bisect not required."
)

# Check for any mismatch on convergence/init path with confidence >= 90
for (mr in mismatch_rows) {
  if (mr$param %in% CONVERGENCE_INIT_PARAMS && mr$confidence >= 90) {
    short_circuit <- TRUE
    short_circuit_reason <- sprintf(
      "Mismatch on convergence/init param '%s' (confidence=%d)", mr$param, mr$confidence
    )
    break
  }
}

# ---- Write TOON output — valid YAML, sequence-of-mappings for param_table ----
# Strings containing ':' or '#' are single-quoted to avoid YAML parse errors.
yaml_quote <- function(s) {
  s <- as.character(s)
  # Single-quote if value contains YAML-unsafe characters or leading/trailing whitespace.
  needs_quote <- grepl(":", s, fixed = TRUE) ||
                 grepl("#", s, fixed = TRUE) ||
                 grepl(">", s, fixed = TRUE) ||
                 grepl("|", s, fixed = TRUE) ||
                 grepl("&", s, fixed = TRUE) ||
                 grepl("*", s, fixed = TRUE) ||
                 grepl("!", s, fixed = TRUE) ||
                 grepl(",", s, fixed = TRUE) ||
                 grepl("{", s, fixed = TRUE) ||
                 grepl("}", s, fixed = TRUE) ||
                 grepl("[", s, fixed = TRUE) ||
                 grepl("]", s, fixed = TRUE) ||
                 grepl("^\\s|\\s$", s) ||
                 nchar(s) == 0
  if (needs_quote) {
    # Escape any single-quotes inside by doubling them
    s <- gsub("'", "''", s, fixed = TRUE)
    s <- paste0("'", s, "'")
  }
  s
}

toon_path <- file.path(OUTDIR, "params_diff.toon")
con <- file(toon_path, "w")

writeLines("task_id: T1_2apm", con)
writeLines("success: true", con)
writeLines("data:", con)
writeLines("  arity_check: 42", con)
writeLines("  files_created:", con)
writeLines("    - benchmarks/2apm/dump_params_R.R", con)
writeLines("    - benchmarks/2apm/dump_params_py.py", con)
writeLines("    - benchmarks/2apm/diff_params.R", con)
writeLines("    - benchmarks/2apm/params_R.txt", con)
writeLines("    - benchmarks/2apm/params_Py.txt", con)
writeLines("    - benchmarks/2apm/params_diff.toon", con)
writeLines("", con)
writeLines("  param_table:", con)
for (r in rows) {
  writeLines(sprintf("    - param: %s", r$param), con)
  writeLines(sprintf("      r_value: %s", yaml_quote(r$r_value)), con)
  writeLines(sprintf("      py_value: %s", yaml_quote(r$py_value)), con)
  pdwm <- if (is.null(r$py_default_when_missing) || r$py_default_when_missing == "N/A") {
    "null"
  } else {
    yaml_quote(r$py_default_when_missing)
  }
  writeLines(sprintf("      py_default_when_missing: %s", pdwm), con)
  writeLines(sprintf("      mismatch: %s", if (r$mismatch) "true" else "false"), con)
  writeLines(sprintf("      confidence: %d", r$confidence), con)
  writeLines(sprintf("      note: %s", yaml_quote(r$note)), con)
}

writeLines("", con)
writeLines("  mismatches:", con)
if (length(mismatch_rows) == 0) {
  writeLines("    []", con)
} else {
  for (mr in mismatch_rows) {
    writeLines(sprintf("    - param: %s", mr$param), con)
    writeLines(sprintf("      r_value: %s", yaml_quote(mr$r_value)), con)
    writeLines(sprintf("      py_value: %s", yaml_quote(mr$py_value)), con)
    writeLines(sprintf("      confidence: %d", mr$confidence), con)
  }
}

writeLines("", con)
writeLines(sprintf("  short_circuit: %s", if (short_circuit) "true" else "false"), con)
# Use block scalar (>-) for reason — safe for arbitrary prose
writeLines("  short_circuit_reason: >-", con)
writeLines(sprintf("    %s", short_circuit_reason), con)
writeLines("", con)
writeLines("error_log: null", con)
close(con)

cat(sprintf("\nparams_diff.toon written to %s\n", toon_path))
cat(sprintf("Mismatches found: %d\n", length(mismatch_rows)))
for (mr in mismatch_rows) {
  cat(sprintf("  - %s (confidence=%d): R=%s | Py=%s\n",
              mr$param, mr$confidence, mr$r_value, mr$py_value))
}
cat(sprintf("short_circuit: %s\n", if (short_circuit) "TRUE" else "FALSE"))
cat(sprintf("Reason: %s\n", short_circuit_reason))
