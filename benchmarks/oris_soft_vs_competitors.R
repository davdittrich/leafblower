# Honest Performance Gate — oris_soft vs. doc-named competitors
# Phase: 03-honest-performance-gate, Plan 01 (tracer slice)
#
# D-04: oris_soft is the single authoritative solver for the headline claim.
# D-07: the competitor set is drawn ONLY from docs/methods/oris.md's
#   "Practitioner implementations & use cases" table — this plan measures
#   `survey::calibrate` (raking), the remaining competitors arrive in plan 02.
# D-08: every number below is measured FRESH in this run. This script does
#   NOT read benchmarks/study/report/tables/*.csv or any other pre-aggregated
#   result file.
# D-05: no comparison against, mention of, or variable named after the
#   unreleased package the old "Nx faster than" framing used as its baseline.
#   That framing is exactly what this phase exists to retire — do not
#   reintroduce it, not even as a commented-out arm.
#
# CSV schema (frozen; plans 02-04 extend ROWS into it, not columns):
#   input_class, n, n_margins, n_categories, m_cell, m_cell_over_n,
#   max_weight, arm, wall_s, max_error, max_w, min_w, deff, n_eff,
#   iterations, ok, note

suppressPackageStartupMessages({
  library(leafblower)
  library(survey)
  library(bench)
})

# --- Determinism guard (CLAUDE.md single-thread BLAS protocol, ENFORCED) ---
# This is the first check the script performs, before any fixture or solve.
# Deliberately does NOT call Sys.setenv() to satisfy itself — the guard only
# has teeth if the caller's environment is what is checked.
require_single_thread_blas <- function() {
  vars <- c("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")
  vals <- Sys.getenv(vars, unset = NA_character_)
  bad  <- vars[is.na(vals) | vals != "1"]
  if (length(bad) > 0L) {
    stop(sprintf(
      "refusing to measure: %s must be set to \"1\" for a reproducible single-thread BLAS run (CLAUDE.md protocol); got {%s}",
      paste(bad, collapse = ", "),
      paste(sprintf("%s=%s", bad, ifelse(is.na(vals[bad]), "<unset>", vals[bad])), collapse = ", ")
    ), call. = FALSE)
  }
}

# Shared, arm-independent accuracy metric — applied identically to every
# arm's returned weight vector. Never reads a package's self-reported
# convergence number. Accumulates the achieved proportion directly (no
# difference-of-large-sums cancellation).
margin_max_error <- function(w, df, tgt) {
  Z <- sum(w)
  max(vapply(names(tgt), function(k) {
    T_k <- tgt[[k]]
    max(vapply(names(T_k), function(j) {
      abs(sum(w[df[[k]] == j]) / Z - T_k[[j]])
    }, numeric(1)))
  }, numeric(1)))
}

# Bound compliance is compared on one scale: sum(w) == n (leafblower's own
# exit convention, lbw::finalize_weights: a single pre-bounds scale to n,
# THEN bounds_mode dispatch). leafblower's returned weights already sit at
# this scale to double-precision (finalize_weights enforces it internally)
# and CLAUDE.md explicitly FORBIDS renormalizing AFTER unit-mode water-fill
# ("silently breaks the bounds_mode='unit' clamps") — confirmed empirically:
# re-scaling leafblower's own output pushed a weight sitting exactly at the
# max_weight bound to 3.0000000036, failing the max_w <= max_weight + 1e-10
# check this script itself asserts. So this helper is applied ONLY to the
# competitor arm below, whose weights (survey's epsilon-tolerance raking)
# do NOT already sit at sum(w) == n, to bring it onto the same scale.
normalize_to_n <- function(w, n) w * n / sum(w)

# One CSV row for one arm. Metrics are passed in pre-computed so each
# call site to margin_max_error() stays visible at the call site, not
# hidden behind this row-shaping helper.
arm_row <- function(input_class, n, n_margins, n_categories, m_cell, max_weight,
                     arm, wall_s, max_error, max_w, min_w, deff, n_eff,
                     iterations, ok, note) {
  data.frame(
    input_class = input_class, n = n, n_margins = n_margins,
    n_categories = n_categories, m_cell = m_cell, m_cell_over_n = m_cell / n,
    max_weight = max_weight, arm = arm, wall_s = wall_s, max_error = max_error,
    max_w = max_w, min_w = min_w, deff = deff, n_eff = n_eff,
    iterations = iterations, ok = ok, note = note,
    stringsAsFactors = FALSE
  )
}

# Runs both arms on one already-built input class (df + tgt) and returns
# the two CSV rows. Plan 02 adds further input classes and competitors by
# calling this same function (or its per-arm blocks) on new fixtures —
# the fixture/measurement split is what makes that extension additive.
run_input_class <- function(input_class, df, tgt, max_weight, n_categories) {
  n <- nrow(df)
  margin_cols <- names(tgt)
  n_margins   <- length(margin_cols)
  # A property of the fixture, identical for every arm — computed once.
  m_cell <- nrow(unique(df[margin_cols]))

  cat(sprintf("\n=== %s n=%d K=%d nj=%d m_cell=%d (m_cell/n=%.4f) ===\n",
              input_class, n, n_margins, n_categories, m_cell, m_cell / n))

  # --- Arm 1: leafblower_oris_soft ---
  # bounds_mode="unit" is load-bearing, not incidental: the default "cell"
  # enforces bounds only in cell aggregate, whereas survey/icarus/ReGenesees
  # all enforce a per-observation bound on the ratio g = w/d. With unit
  # design weights, bounds_mode="unit" + max_weight=3 corresponds exactly to
  # the competitors' bounds=c(0,3) on g — comparing a cell-aggregate
  # leafblower arm against per-observation competitors would be a different
  # optimisation problem and the timing comparison would be meaningless.
  # convergence is oris_soft's own CANONICAL default (R/harvest.R roxygen,
  # `convergence` param): metric="marginal_kl", rule="improvement", tol=0.001.
  # Passed explicitly (not convergence=list()) so the intent is visible at
  # the call site. Originally this arm passed the `absolute` shorthand,
  # which forces metric="max_err"/rule="threshold" — survey::calibrate's
  # own native stopping rule, not oris_soft's. That silently made oris_soft
  # stop the instant it first crossed the tolerance instead of continuing
  # to refine on its own plateau-detection rule, exactly the "stopped
  # early" confound this phase's determinism protocol exists to eliminate.
  # Fixed per tracer-checkpoint review; see 03-01-SUMMARY.md.
  lb_call <- function() {
    harvest(df, tgt, method = "oris_soft", max_weight = max_weight,
            bounds_mode = "unit", attach_weights = FALSE,
            convergence = list(metric = "marginal_kl", rule = "improvement", tol = 0.001))
  }
  bm_lb  <- bench::mark(run = lb_call(), iterations = 2, check = FALSE,
                         memory = FALSE, filter_gc = FALSE)
  w_lb   <- lb_call()
  res_lb <- attr(w_lb, "result")
  # NOT normalized — see the comment on normalize_to_n() above.
  w_lb_n <- as.numeric(w_lb)
  max_error_lb <- margin_max_error(w_lb_n, df, tgt)
  # status 0 = converged; 5 = plateau at constrained optimum (weights valid,
  # the bound is legitimately active) — both are a usable result. 1-4 are not.
  ok_lb <- isTRUE(res_lb$status %in% c(0L, 5L))
  note_lb <- sprintf("convergence=list(metric='marginal_kl',rule='improvement',tol=0.001) (oris_soft canonical default) requested; status=%d, iterations=%d",
                      res_lb$status, res_lb$iterations)
  row_lb <- arm_row(input_class, n, n_margins, n_categories, m_cell, max_weight,
                     "leafblower_oris_soft", as.numeric(bm_lb$median), max_error_lb,
                     max(w_lb_n), min(w_lb_n),
                     leafblower::design_effect(w_lb_n),
                     leafblower::effective_sample_size(w_lb_n),
                     res_lb$iterations, ok_lb, note_lb)
  cat(sprintf("  %-22s wall=%7.4fs status=%d max_err=%.3e max_w=%.3f n_eff=%.1f\n",
              "leafblower_oris_soft", row_lb$wall_s, res_lb$status,
              row_lb$max_error, row_lb$max_w, row_lb$n_eff))

  # --- Arm 2: survey_calibrate ---
  # Population totals derived from the SAME tgt list, scaled to n, so both
  # arms are given identical control totals (confound control). formula/
  # population as lists-of-formulas/lists-of-tables — the classical-raking
  # special case documented under ?survey::calibrate ("in the same format
  # as the input to rake"); read from the installed help, not assumed.
  formula_list <- lapply(margin_cols, function(k) stats::as.formula(paste0("~", k)))
  population_list <- lapply(margin_cols, function(k) {
    T_k <- tgt[[k]]
    data.frame(setNames(list(names(T_k)), k), Freq = as.numeric(T_k) * n)
  })
  design <- survey::svydesign(ids = ~1, weights = ~1, data = df)
  sv_call <- function() {
    survey::calibrate(design, formula_list, population_list, calfun = "raking",
                       bounds = c(0, max_weight), epsilon = 1e-3)
  }
  # A single competitor failure must not cost the maintainer the whole run
  # (task 3): every competitor arm's call is tryCatch()-isolated, a failed
  # arm emits ok=FALSE with the condition message in `note`, and the other
  # arms still run.
  row_sv <- tryCatch({
    bm_sv  <- bench::mark(run = sv_call(), iterations = 2, check = FALSE,
                           memory = FALSE, filter_gc = FALSE)
    cal_sv <- sv_call()
    w_sv   <- stats::weights(cal_sv)
    w_sv_n <- normalize_to_n(as.numeric(w_sv), n)
    max_error_sv <- margin_max_error(w_sv_n, df, tgt)
    ok_sv <- all(is.finite(w_sv)) && max(w_sv) <= max_weight + 1e-6
    note_sv <- "epsilon=1e-3 requested (survey::calibrate, calfun='raking', bounds=c(0,max_weight))"
    # iterations is not comparable across packages for this arm.
    arm_row(input_class, n, n_margins, n_categories, m_cell, max_weight,
            "survey_calibrate", as.numeric(bm_sv$median), max_error_sv,
            max(w_sv_n), min(w_sv_n),
            leafblower::design_effect(w_sv_n),
            leafblower::effective_sample_size(w_sv_n),
            NA_integer_, ok_sv, note_sv)
  }, error = function(e) {
    arm_row(input_class, n, n_margins, n_categories, m_cell, max_weight,
            "survey_calibrate", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
            NA_integer_, FALSE, paste0("error: ", substr(conditionMessage(e), 1, 200)))
  })
  cat(sprintf("  %-22s wall=%7s max_err=%9s max_w=%7s n_eff=%9s\n",
              "survey_calibrate",
              if (is.finite(row_sv$wall_s)) sprintf("%.4fs", row_sv$wall_s) else "NA",
              if (is.finite(row_sv$max_error)) sprintf("%.3e", row_sv$max_error) else "NA",
              if (is.finite(row_sv$max_w)) sprintf("%.3f", row_sv$max_w) else "NA",
              if (is.finite(row_sv$n_eff)) sprintf("%.1f", row_sv$n_eff) else "NA"))

  # --- Arm 3: icarus_calibration ---
  # D-09 dependency scoping: requireNamespace() guard, not a DESCRIPTION
  # Suggests: entry. icarus's bounds argument is only honoured by its
  # "logit" method (per ?icarus::calibration: "bounds: ... for bounded
  # methods ('logit')") -- logit is the bounded-calibration analog to
  # survey's calfun="raking"/bounds=c(0,max_weight), read from the
  # installed help before writing this call. Margin totals are entered via
  # newMarginMatrix()/addMargin() (magrittr-free form) in the SAME level
  # order as tgt[[k]] (setNames(p_skew, levels(df[[k]])) in the medium
  # fixture below), which is the order icarus::calibrationMatrix()
  # internally assigns to dummy columns (verified: data.matrix() on a
  # factor -> ascending integer codes -> colToDummies() in that order).
  if (!requireNamespace("icarus", quietly = TRUE)) {
    row_icarus <- arm_row(input_class, n, n_margins, n_categories, m_cell, max_weight,
                           "icarus_calibration", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                           NA_integer_, FALSE,
                           "skipped: requireNamespace('icarus') returned FALSE; install with install.packages('icarus') to include this arm")
  } else {
    row_icarus <- tryCatch({
      df_ic <- df
      df_ic$w0_icarus <- 1
      mm_ic <- icarus::newMarginMatrix()
      for (k in margin_cols) mm_ic <- icarus::addMargin(mm_ic, k, as.numeric(tgt[[k]]))
      ic_call <- function() {
        icarus::calibration(data = df_ic, marginMatrix = mm_ic, colWeights = "w0_icarus",
                             method = "logit", bounds = c(0, max_weight), popTotal = n,
                             pct = TRUE, description = FALSE, calibTolerance = 1e-6,
                             maxIter = 2500)
      }
      bm_ic <- bench::mark(run = ic_call(), iterations = 2, check = FALSE,
                            memory = FALSE, filter_gc = FALSE)
      w_ic   <- ic_call()
      w_ic_n <- normalize_to_n(as.numeric(w_ic), n)
      max_error_ic <- margin_max_error(w_ic_n, df, tgt)
      ok_ic <- all(is.finite(w_ic_n)) && max(w_ic_n) <= max_weight + 1e-6
      note_ic <- sprintf(
        "method='logit', bounds=c(0,%g), calibTolerance=1e-6 requested (icarus::calibration; logit is icarus's bounds-supporting method per its own help)",
        max_weight)
      arm_row(input_class, n, n_margins, n_categories, m_cell, max_weight,
              "icarus_calibration", as.numeric(bm_ic$median), max_error_ic,
              max(w_ic_n), min(w_ic_n),
              leafblower::design_effect(w_ic_n),
              leafblower::effective_sample_size(w_ic_n),
              NA_integer_, ok_ic, note_ic)
    }, error = function(e) {
      arm_row(input_class, n, n_margins, n_categories, m_cell, max_weight,
              "icarus_calibration", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
              NA_integer_, FALSE, paste0("error: ", substr(conditionMessage(e), 1, 200)))
    })
  }
  cat(sprintf("  %-22s wall=%7s max_err=%9s max_w=%7s n_eff=%9s\n",
              "icarus_calibration",
              if (is.finite(row_icarus$wall_s)) sprintf("%.4fs", row_icarus$wall_s) else "NA",
              if (is.finite(row_icarus$max_error)) sprintf("%.3e", row_icarus$max_error) else "NA",
              if (is.finite(row_icarus$max_w)) sprintf("%.3f", row_icarus$max_w) else "NA",
              if (is.finite(row_icarus$n_eff)) sprintf("%.1f", row_icarus$n_eff) else "NA"))

  # --- Arm 4: ReGenesees_e_calibrate ---
  # D-09 dependency scoping: requireNamespace() guard, not a DESCRIPTION
  # Suggests: entry. Population totals go through ReGenesees's own
  # pop.template()/fill.template() contract (read from the installed help
  # before writing this call): pop.template() builds a totals-slot data
  # frame from a calmodel formula; this script has no sampling-frame
  # microdata to hand fill.template() (only aggregate proportions in tgt),
  # so the template's NA slots are filled directly by name -- template
  # column names are "<margin><level>" (e.g. "m1a"); parsed back into
  # margin/level and looked up in tgt. calmodel is additive
  # ("~m1+m2+...-1", the classical-raking form for independently-specified
  # margins), matching survey/icarus's per-margin (not joint) targets.
  if (!requireNamespace("ReGenesees", quietly = TRUE)) {
    row_regenesees <- arm_row(input_class, n, n_margins, n_categories, m_cell, max_weight,
                               "ReGenesees_e_calibrate", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                               NA_integer_, FALSE,
                               "skipped: requireNamespace('ReGenesees') returned FALSE; install with install.packages('ReGenesees') to include this arm")
  } else {
    row_regenesees <- tryCatch({
      df_rg <- df
      df_rg$w0_regenesees <- 1
      df_rg$id_regenesees <- seq_len(n)
      calmodel_rg <- stats::as.formula(paste("~", paste(margin_cols, collapse = "+"), "-1"))
      tmpl_rg <- ReGenesees::pop.template(data = df_rg, calmodel = calmodel_rg)
      margin_cols_by_len <- margin_cols[order(-nchar(margin_cols))]
      for (cn in names(tmpl_rg)) {
        mcol <- margin_cols_by_len[startsWith(cn, margin_cols_by_len)][1]
        lvl  <- substring(cn, nchar(mcol) + 1)
        tmpl_rg[1, cn] <- tgt[[mcol]][[lvl]] * n
      }
      rg_call <- function() {
        des_rg <- ReGenesees::e.svydesign(data = df_rg, ids = ~id_regenesees,
                                           weights = ~w0_regenesees)
        cal_rg <- ReGenesees::e.calibrate(des_rg, tmpl_rg, calmodel = calmodel_rg,
                                           calfun = "raking", bounds = c(0, max_weight))
        stats::weights(cal_rg)
      }
      bm_rg <- bench::mark(run = rg_call(), iterations = 2, check = FALSE,
                            memory = FALSE, filter_gc = FALSE)
      w_rg   <- rg_call()
      w_rg_n <- normalize_to_n(as.numeric(w_rg), n)
      max_error_rg <- margin_max_error(w_rg_n, df, tgt)
      ok_rg <- all(is.finite(w_rg_n)) && max(w_rg_n) <= max_weight + 1e-6
      note_rg <- sprintf(
        "calmodel=%s, calfun='raking', bounds=c(0,%g) requested (ReGenesees::e.calibrate; pop.template()/manual fill, no sampling-frame fill.template() available)",
        deparse(calmodel_rg), max_weight)
      arm_row(input_class, n, n_margins, n_categories, m_cell, max_weight,
              "ReGenesees_e_calibrate", as.numeric(bm_rg$median), max_error_rg,
              max(w_rg_n), min(w_rg_n),
              leafblower::design_effect(w_rg_n),
              leafblower::effective_sample_size(w_rg_n),
              NA_integer_, ok_rg, note_rg)
    }, error = function(e) {
      arm_row(input_class, n, n_margins, n_categories, m_cell, max_weight,
              "ReGenesees_e_calibrate", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
              NA_integer_, FALSE, paste0("error: ", substr(conditionMessage(e), 1, 200)))
    })
  }
  cat(sprintf("  %-22s wall=%7s max_err=%9s max_w=%7s n_eff=%9s\n",
              "ReGenesees_e_calibrate",
              if (is.finite(row_regenesees$wall_s)) sprintf("%.4fs", row_regenesees$wall_s) else "NA",
              if (is.finite(row_regenesees$max_error)) sprintf("%.3e", row_regenesees$max_error) else "NA",
              if (is.finite(row_regenesees$max_w)) sprintf("%.3f", row_regenesees$max_w) else "NA",
              if (is.finite(row_regenesees$n_eff)) sprintf("%.1f", row_regenesees$n_eff) else "NA"))

  rbind(row_lb, row_sv, row_icarus, row_regenesees)
}

# --- Fixture: medium_100k_5margins ---
# nj^K = 4^5 = 1024, n = 100000 -> m_cell/n far below 1: a non-degenerate,
# compression-benefiting class, deliberately NOT the K=20 uniform-random
# class RESEARCH.md Pitfall 3 warns is the wrong fixture for a positive claim.
# This shape matches the PRD's own contradictory medium-scale target
# (100K rows, 5 margins, categories per margin) so this number directly
# retires that contradiction (SC2) rather than sitting beside it.
require_single_thread_blas()

set.seed(304)
n  <- 100000L
K  <- 5L
nj <- 4L
margin_cols <- paste0("m", seq_len(K))
df <- as.data.frame(lapply(seq_len(K), function(k)
  factor(sample(letters[seq_len(nj)], n, TRUE))))
names(df) <- margin_cols

# Skewed, non-uniform per-margin target so calibration is a real adjustment,
# not a no-op. NOTE: the neighbourhood of 0.45/0.30/0.15/0.10 named in the
# plan is infeasible at max_weight=3 for K=5 INDEPENDENT margins — the
# combined multiplicative correction needed on a cell disfavoured across
# all 5 margins exceeds the bound and both arms fail to converge
# (oris_soft: status=4 budget-exhausted at max_error=3.8e-2; survey:
# "Calibration failed"). 0.40/0.28/0.18/0.14 is also too tight at n=100000
# (oris_soft hits status=5 constrained-optimum plateau at max_error=1.49e-3,
# above the 1e-3 floor task 3 gates on). 0.36/0.27/0.20/0.17 keeps the joint
# correction inside the max_weight=3 bound (leafblower still clamps weights
# at the bound on this fixture — the bound is exercised, not a no-op) while
# both arms converge (status=0) well under the 1e-3 accuracy floor. See
# 03-01-SUMMARY.md for the measured evidence at each candidate skew.
p_skew <- c(0.36, 0.27, 0.20, 0.17)
p_skew <- p_skew / sum(p_skew)
tgt <- setNames(lapply(margin_cols, function(k) setNames(p_skew, levels(df[[k]]))), margin_cols)

results <- run_input_class("medium_100k_5margins", df, tgt, max_weight = 3, n_categories = nj)

# --- G-03-1 gap closure: leafblower-only arms for oris/raking/newton_kl ---
# oris/raking/newton_kl all attack the identical K-margin box-bounded KL
# calibration problem survey_calibrate/icarus_calibration/ReGenesees_e_calibrate
# already measured above on this same fixture (docs/methods/oris.md,
# docs/methods/raking.md, docs/methods/newton_kl.md's own Practitioner-
# implementations tables all name survey::calibrate as the direct competitor
# for this objective class) -- D-08's "narrow, not a duplicate study" framing
# means those three competitor rows are reused here, not recomputed a second
# time for each new leafblower method.
m_cell_medium <- nrow(unique(df[margin_cols]))
lb_only_arm_row <- function(method_name, arm_label) {
  lb_call <- function() {
    harvest(df, tgt, method = method_name, max_weight = 3,
            bounds_mode = "unit", attach_weights = FALSE, convergence = list())
  }
  bm_lb  <- bench::mark(run = lb_call(), iterations = 2, check = FALSE,
                         memory = FALSE, filter_gc = FALSE)
  w_lb   <- lb_call()
  res_lb <- attr(w_lb, "result")
  # NOT normalized -- leafblower's own arms are never renormalized (see the
  # comment on normalize_to_n() above).
  w_lb_n <- as.numeric(w_lb)
  max_error_lb <- margin_max_error(w_lb_n, df, tgt)
  ok_lb <- isTRUE(res_lb$status %in% c(0L, 5L))
  note_lb <- sprintf(
    "convergence=list() (per-method natural default per R/harvest.R:424); status=%d, iterations=%d",
    res_lb$status, res_lb$iterations)
  row <- arm_row("medium_100k_5margins", n, K, nj, m_cell_medium,
                  3, arm_label, as.numeric(bm_lb$median), max_error_lb,
                  max(w_lb_n), min(w_lb_n),
                  leafblower::design_effect(w_lb_n),
                  leafblower::effective_sample_size(w_lb_n),
                  res_lb$iterations, ok_lb, note_lb)
  cat(sprintf("  %-22s wall=%7.4fs status=%d max_err=%.3e max_w=%.3f n_eff=%.1f\n",
              arm_label, row$wall_s, res_lb$status, row$max_error, row$max_w, row$n_eff))
  row
}

results <- rbind(results,
                  lb_only_arm_row("oris", "leafblower_oris"),
                  lb_only_arm_row("raking", "leafblower_raking"),
                  lb_only_arm_row("newton_kl", "leafblower_newton_kl"))

# --- G-03-1 gap closure: leafblower_greg / leafblower_logit vs distance-
# matched NEW survey::calibrate(calfun=) competitor rows ---
# greg (chi-square/linear distance, docs/methods/greg.md) and logit (logit
# distance, docs/methods/logit.md, which names survey::calibrate(calfun=
# "logit") first in its own Practitioner-implementations table) attack
# DIFFERENT objectives than the raking-calfun arm already in `results` --
# pairing them against that arm would be an objective mismatch, so each
# gets its own calfun-matched survey::calibrate row. Reuses the medium
# fixture's already-built design/formula_list/population_list (from the
# earlier survey_calibrate block inside run_input_class(), rebuilt
# identically here since those objects are local to that function).
formula_list_medium <- lapply(margin_cols, function(k) stats::as.formula(paste0("~", k)))
population_list_medium <- lapply(margin_cols, function(k) {
  T_k <- tgt[[k]]
  data.frame(setNames(list(names(T_k)), k), Freq = as.numeric(T_k) * n)
})
design_medium <- survey::svydesign(ids = ~1, weights = ~1, data = df)

survey_calfun_arm_row <- function(calfun_name, arm_label) {
  tryCatch({
    sv_call <- function() {
      survey::calibrate(design_medium, formula_list_medium, population_list_medium,
                         calfun = calfun_name, bounds = c(0, 3))
    }
    bm_sv  <- bench::mark(run = sv_call(), iterations = 2, check = FALSE,
                           memory = FALSE, filter_gc = FALSE)
    cal_sv <- sv_call()
    w_sv   <- stats::weights(cal_sv)
    w_sv_n <- normalize_to_n(as.numeric(w_sv), n)
    max_error_sv <- margin_max_error(w_sv_n, df, tgt)
    ok_sv <- all(is.finite(w_sv)) && max(w_sv) <= 3 + 1e-6
    note_sv <- sprintf("calfun='%s' requested (survey::calibrate, bounds=c(0,3))", calfun_name)
    row <- arm_row("medium_100k_5margins", n, K, nj, m_cell_medium, 3,
                    arm_label, as.numeric(bm_sv$median), max_error_sv,
                    max(w_sv_n), min(w_sv_n),
                    leafblower::design_effect(w_sv_n),
                    leafblower::effective_sample_size(w_sv_n),
                    NA_integer_, ok_sv, note_sv)
    cat(sprintf("  %-22s wall=%7.4fs max_err=%.3e max_w=%.3f n_eff=%.1f\n",
                arm_label, row$wall_s, row$max_error, row$max_w, row$n_eff))
    row
  }, error = function(e) {
    arm_row("medium_100k_5margins", n, K, nj, m_cell_medium, 3,
            arm_label, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
            NA_integer_, FALSE, paste0("error: ", substr(conditionMessage(e), 1, 200)))
  })
}

row_lb_greg <- lb_only_arm_row("greg", "leafblower_greg")
row_lb_logit <- lb_only_arm_row("logit", "leafblower_logit")
# icarus_calibration (already computed on this fixture, internally
# method='logit' per its own D-09-scoped note above) is ALSO a same-distance
# match for this arm -- pointing to its existing row, not recomputed.
row_lb_logit$note <- paste0(
  row_lb_logit$note,
  "; icarus_calibration (see its own row on this fixture) is also a same-distance (logit) match for this arm")

results <- rbind(results, row_lb_greg, row_lb_logit,
                  survey_calfun_arm_row("linear", "survey_calibrate_linear"),
                  survey_calfun_arm_row("logit", "survey_calibrate_logit"))

# --- G-03-1 gap closure: leafblower_chebyshev vs optweight::optweight.svy(norm='linf') ---
# chebyshev is the last R-side method needing a doc-named competitor not
# already installed for this phase (docs/methods/chebyshev.md's own
# Practitioner-implementations table: "The optweight package is the primary
# open-source tool supporting the minimax norm explicitly"). Package
# legitimacy verified this session (CRAN PDF manual fetch, 03-RESEARCH.md's
# extended Package Legitimacy Audit): optweight v2.0.1, Noah Greifer (also
# maintains WeightIt/cobalt), actively maintained, genuinely supports
# norm="linf" (f(w,b,s) = max_i s_i|w_i-b_i|).
#
# convergence=list() -- chebyshev ignores rule/metric/stop_when entirely
# (R/harvest.R's own documented parameter note: its interior-point solver
# stops on its own complementarity-gap criterion), so lb_only_arm_row()'s
# existing convergence=list() default applies unmodified -- no new helper
# needed for this arm.
row_lb_chebyshev <- lb_only_arm_row("chebyshev", "leafblower_chebyshev")

# optweight_linf: optweight.svy() is the "Stable Balancing Weights for
# Generalization" variant (no treat= argument) -- the shape matching a pure
# population-calibration problem with no binary/multi-category treatment
# variable, unlike plain optweight() whose primary examples are all
# treatment-balancing. targets built via unlist(tgt[margin_cols],
# use.names=FALSE): empirically verified against the installed package this
# session (process_targets() run on this exact fixture shape) to require
# ALL levels per factor margin in level order, NOT non-reference-category-
# only -- process_targets()'s own Rd ("For factor variables, a target value
# must be specified for each level of the factor, and these values must add
# up to 1") and optweight's internal model.covs expansion (m1_a,m1_b,m1_c,
# m1_d,m2_a,...) confirm every level is a separate dummy column, none
# dropped as a reference. tols=0 requests exact margin matching (closest
# analogue to harvest()'s own equality-constraint margins). min.w=1e-8
# mirrors the package's own documented default.
if (!requireNamespace("optweight", quietly = TRUE)) {
  row_ow <- arm_row("medium_100k_5margins", n, K, nj, m_cell_medium, 3,
                     "optweight_linf", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                     NA_integer_, FALSE,
                     "skipped: requireNamespace('optweight') returned FALSE; install with install.packages('optweight') to include this arm")
} else {
  row_ow <- tryCatch({
    formula_ow  <- stats::reformulate(margin_cols)
    targets_ow  <- unlist(tgt[margin_cols], use.names = FALSE)
    ow_call <- function() {
      optweight::optweight.svy(formula_ow, data = df, targets = targets_ow,
                                tols = 0, norm = "linf", min.w = 1e-8)
    }
    bm_ow  <- bench::mark(run = ow_call(), iterations = 2, check = FALSE,
                           memory = FALSE, filter_gc = FALSE)
    cal_ow <- ow_call()
    w_ow   <- stats::weights(cal_ow)
    w_ow_n <- normalize_to_n(as.numeric(w_ow), n)
    max_error_ow <- margin_max_error(w_ow_n, df, tgt)
    note_ow <- paste0(
      "optweight.svy(norm='linf') minimizes max_i s_i|w_i-b_i| (deviation from base weight, default 1), ",
      "margin balance enforced via tols=0 as a CONSTRAINT, not the objective; leafblower's chebyshev instead ",
      "minimizes max weighted margin error DIRECTLY as the LP objective (weights unconstrained in the objective) -- ",
      "a related but non-identical minimax formulation, per docs/methods/chebyshev.md's own citation of optweight ",
      "as the primary open-source tool supporting the minimax norm explicitly. optweight has no max.w argument ",
      "(only min.w, confirmed against the CRAN manual fetched 2026-08-15) so max_w bound-compliance is ",
      "UNVERIFIABLE on this arm.")
    # No max.w argument exists on this package (confirmed against
    # optweight.svy's own formals: formula, data, tols, targets, s.weights,
    # b.weights, norm, min.w, verbose, ...) -- ok is NA (not TRUE/FALSE),
    # not claiming a bound compliance the package cannot express; max_w
    # reports the raw achieved maximum, on the same normalize_to_n()
    # n-scale every other competitor row in this file uses for
    # comparability, not as a compliance claim.
    arm_row("medium_100k_5margins", n, K, nj, m_cell_medium, 3,
            "optweight_linf", as.numeric(bm_ow$median), max_error_ow,
            max(w_ow_n), min(w_ow_n),
            leafblower::design_effect(w_ow_n),
            leafblower::effective_sample_size(w_ow_n),
            NA_integer_, NA, note_ow)
  }, error = function(e) {
    arm_row("medium_100k_5margins", n, K, nj, m_cell_medium, 3,
            "optweight_linf", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
            NA_integer_, FALSE, paste0("error: ", substr(conditionMessage(e), 1, 200)))
  })
}
cat(sprintf("  %-22s wall=%7s max_err=%9s max_w=%7s n_eff=%9s\n",
            "optweight_linf",
            if (is.finite(row_ow$wall_s)) sprintf("%.4fs", row_ow$wall_s) else "NA",
            if (is.finite(row_ow$max_error)) sprintf("%.3e", row_ow$max_error) else "NA",
            if (is.finite(row_ow$max_w)) sprintf("%.3f", row_ow$max_w) else "NA",
            if (is.finite(row_ow$n_eff)) sprintf("%.1f", row_ow$n_eff) else "NA"))

results <- rbind(results, row_lb_chebyshev, row_ow)

# --- Fixture: large_stepstone_fulldata ---
# The tracked 1,582,732-row / 9-margin / 836-category real-survey fixture
# (benchmarks/stepstone_fulldata_bench_data.parquet / _targets.json), reused
# per D-08's "reuse benchmarks/ infrastructure, do not build a parallel
# harness" instruction — this is stepstone_fulldata_benchmark.R's own
# fixture, not a newly-invented large fixture. SC1's large-scale figure
# lives here. Loading follows that script's convention (arrow::read_parquet
# + jsonlite::fromJSON); its comparison target (autumn) and its
# MAX_WEIGHT=5/method="oris"/default bounds_mode are NOT reused — this arm
# uses max_weight=3/bounds_mode="unit" to match the medium class's bound
# convention, so both leafblower figures can be quoted in the same sentence.
large_parquet_path <- "benchmarks/stepstone_fulldata_bench_data.parquet"
large_targets_path <- "benchmarks/stepstone_fulldata_bench_targets.json"

if (!file.exists(large_parquet_path) || !file.exists(large_targets_path)) {
  missing_large_file <- if (!file.exists(large_parquet_path)) large_parquet_path else large_targets_path
  cat(sprintf("\n=== large_stepstone_fulldata: SKIPPED (missing %s) ===\n", missing_large_file))
  results <- rbind(results, arm_row(
    "large_stepstone_fulldata", NA_integer_, NA_integer_, NA_integer_, NA_integer_, 3,
    "leafblower_oris_soft", NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
    NA_integer_, FALSE, sprintf("skipped: fixture file not found: %s", missing_large_file)))
} else {
  df_large  <- as.data.frame(arrow::read_parquet(large_parquet_path))
  # JSON round-tripping a per-margin proportion vector drifts a few units in
  # the last few decimal places (sums observed at 0.9994-1.0001, not exactly
  # 1 within harvest()'s 1e-6 tolerance) -- renormalise on load, the same
  # correction stepstone_fulldata_benchmark.R itself applies after dropping
  # missing cells (`target_anes[[nm]] <- tgt / sum(tgt)`).
  tgt_large <- lapply(jsonlite::fromJSON(large_targets_path), function(x) {
    v <- unlist(x)
    v / sum(v)
  })
  # Select margin columns by name intersection, not position: a fixture
  # regeneration that reorders columns cannot silently mis-map a margin to
  # the wrong target.
  margin_cols_large <- intersect(names(df_large), names(tgt_large))
  n_large <- nrow(df_large)
  n_margins_large <- length(margin_cols_large)
  # Heterogeneous category counts per margin (2..408 here) -- n_categories
  # is the TOTAL across margins for this class, unlike the uniform
  # medium/known-limit classes where every margin shares one category count.
  n_categories_large <- sum(vapply(tgt_large[margin_cols_large], length, integer(1)))
  m_cell_large <- nrow(unique(df_large[margin_cols_large]))

  cat(sprintf("\n=== large_stepstone_fulldata n=%d K=%d n_categories=%d m_cell=%d (m_cell/n=%.4f) ===\n",
              n_large, n_margins_large, n_categories_large, m_cell_large, m_cell_large / n_large))
  cat("  leafblower_oris_soft   solving...\n")

  lb_large_call <- function() {
    harvest(df_large, tgt_large[margin_cols_large], method = "oris_soft", max_weight = 3,
            bounds_mode = "unit", attach_weights = FALSE,
            convergence = list(metric = "marginal_kl", rule = "improvement", tol = 0.001))
  }
  bm_lb_large  <- bench::mark(run = lb_large_call(), iterations = 2, check = FALSE,
                               memory = FALSE, filter_gc = FALSE)
  w_lb_large   <- lb_large_call()
  res_lb_large <- attr(w_lb_large, "result")
  w_lb_large_n <- as.numeric(w_lb_large)
  max_error_lb_large <- margin_max_error(w_lb_large_n, df_large, tgt_large[margin_cols_large])
  # status 0 = converged; 5 = plateau at constrained optimum -- both usable, as
  # in the medium class (many cells here are legitimately water-filled to the
  # bound given 688 flagged sparse categories in this fixture).
  ok_lb_large <- isTRUE(res_lb_large$status %in% c(0L, 5L))
  note_lb_large <- sprintf(
    "convergence=list(metric='marginal_kl',rule='improvement',tol=0.001) (oris_soft canonical default) requested; status=%d, iterations=%d",
    res_lb_large$status, res_lb_large$iterations)
  row_lb_large <- arm_row("large_stepstone_fulldata", n_large, n_margins_large,
                           n_categories_large, m_cell_large, 3,
                           "leafblower_oris_soft", as.numeric(bm_lb_large$median),
                           max_error_lb_large, max(w_lb_large_n), min(w_lb_large_n),
                           leafblower::design_effect(w_lb_large_n),
                           leafblower::effective_sample_size(w_lb_large_n),
                           res_lb_large$iterations, ok_lb_large, note_lb_large)
  cat(sprintf("  %-22s wall=%7.4fs status=%d max_err=%.3e max_w=%.3f n_eff=%.1f\n",
              "leafblower_oris_soft", row_lb_large$wall_s, res_lb_large$status,
              row_lb_large$max_error, row_lb_large$max_w, row_lb_large$n_eff))

  # Competitors are deliberately NOT run at this scale. survey::calibrate,
  # icarus::calibration and ReGenesees::e.calibrate each build one dense
  # observation-by-category model matrix; at this n and n_categories that
  # matrix alone projects to n * n_categories * 8 bytes. Computed here from
  # the actual fixture shape, not asserted -- an honest, checkable
  # feasibility boundary and itself a publishable comparative finding, not a
  # silent omission of the competitors from the large-scale claim.
  matrix_bytes_large <- as.numeric(n_large) * as.numeric(n_categories_large) * 8
  matrix_gb_large <- matrix_bytes_large / 1024^3
  competitor_note_large <- sprintf(
    "skipped: dense obs-by-category model matrix at n=%d x n_categories=%d projects to ~%.1f GB (n*n_categories*8 bytes); infeasible at this scale",
    n_large, n_categories_large, matrix_gb_large)
  competitor_rows_large <- do.call(rbind, lapply(
    c("survey_calibrate", "icarus_calibration", "ReGenesees_e_calibrate"),
    function(arm_name) arm_row("large_stepstone_fulldata", n_large, n_margins_large,
                                n_categories_large, m_cell_large, 3, arm_name,
                                NA_real_, NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                                NA_integer_, FALSE, competitor_note_large)))
  cat(sprintf("  competitors skipped: %s\n", competitor_note_large))

  results <- rbind(results, row_lb_large, competitor_rows_large)
}

# --- Fixture: known_limit_k20_uniform ---
# Fixture PARAMETERS (n, K, category count per margin, seed, max_weight,
# column naming m1..m20) are byte-identical to
# tests/testthat/test-bench-gate.R's kk1204 block: set.seed(1204),
# n=500000, K=20 uniform-random categorical columns each drawn from 5
# levels, max_weight=3. This class exists to make SC3's known limit a
# MEASURED m_cell_over_n rather than a citation: at K=20 independent
# 5-level columns and n=500,000, effectively every row is its own cell
# (measured below), so ORIS-family cell-compression yields zero benefit.
#
# TARGET distribution deliberately does NOT reuse that test block's
# uniform rep(1/cats, cats): measured on this exact data (see
# 03-02-SUMMARY.md), uniform targets on uniformly-sampled data are
# trivially fittable (max_error ~4e-15, <10 iterations) -- they demonstrate
# zero-compression-benefit but NOT any accuracy ceiling, so they cannot
# back D-01/SC3's "known limit is unachievable" claim, which is about
# accuracy under bounds, not raw cell compression. This class instead
# reuses the ORIGINAL kk1204 investigation's skewed target (0.3, 0.175,
# 0.175, 0.175, 0.175 per margin -- docs/investigations/2026-04-23-kk1204-
# convergence.md "Input" section), the actual parameterization that
# produced the near-infeasible plateau D-01 cites. This is a discrepancy
# from the in-repo test's literal target values, reported here per this
# task's own read_first instruction ("the measured number this task
# produces must be consistent with [the investigation] or the discrepancy
# must be reported") rather than silently forced to match a fixture that
# would not demonstrate the finding this class exists to measure.
set.seed(1204)
n_kk    <- 500000L
K_kk    <- 20L
cats_kk <- 5L
margin_cols_kk <- paste0("m", seq_len(K_kk))
df_kk <- as.data.frame(lapply(seq_len(K_kk), function(k)
  factor(sample(seq_len(cats_kk), n_kk, replace = TRUE))))
names(df_kk) <- margin_cols_kk
skew_kk <- c(0.3, 0.175, 0.175, 0.175, 0.175)
tgt_kk  <- setNames(lapply(margin_cols_kk, function(k)
  setNames(skew_kk, as.character(seq_len(cats_kk)))), margin_cols_kk)

m_cell_kk <- nrow(unique(df_kk[margin_cols_kk]))
cat(sprintf("\n=== known_limit_k20_uniform n=%d K=%d m_cell=%d (m_cell/n=%.4f) ===\n",
            n_kk, K_kk, m_cell_kk, m_cell_kk / n_kk))

# Bounded cost: the investigation's extrapolation is ~1.76 million
# iterations to reach the retired 1e-6 target -- roughly two weeks at this
# host's measured per-iteration cost. A budget in the same order as the
# existing kk1204 test block's 500 caps the run; the achieved max_error at
# that budget IS the finding, not a target this class chases.
kk_iter_budget <- 500L

cat("  leafblower_oris_soft          solving...\n")
kk_lb_call <- function() {
  harvest(df_kk, tgt_kk, method = "oris_soft", max_weight = 3, bounds_mode = "unit",
          attach_weights = FALSE, max_iterations = kk_iter_budget,
          convergence = list(metric = "marginal_kl", rule = "improvement", tol = 0.001))
}
bm_kk_lb  <- bench::mark(run = kk_lb_call(), iterations = 2, check = FALSE,
                          memory = FALSE, filter_gc = FALSE)
w_kk_lb   <- kk_lb_call()
res_kk_lb <- attr(w_kk_lb, "result")
w_kk_lb_n <- as.numeric(w_kk_lb)
max_error_kk_lb <- margin_max_error(w_kk_lb_n, df_kk, tgt_kk)
ok_kk_lb <- isTRUE(res_kk_lb$status %in% c(0L, 5L))
note_kk_lb <- sprintf(
  "max_iterations=%d budget (kk1204-order-of-magnitude cap per docs/investigations/2026-04-23-kk1204-convergence.md); status=%d, iterations=%d; skewed target (0.3/0.175x4 per that investigation), NOT test-bench-gate.R's uniform 1/5 -- see 03-02-SUMMARY.md",
  kk_iter_budget, res_kk_lb$status, res_kk_lb$iterations)
row_kk_lb <- arm_row("known_limit_k20_uniform", n_kk, K_kk, cats_kk, m_cell_kk, 3,
                      "leafblower_oris_soft", as.numeric(bm_kk_lb$median), max_error_kk_lb,
                      max(w_kk_lb_n), min(w_kk_lb_n),
                      leafblower::design_effect(w_kk_lb_n),
                      leafblower::effective_sample_size(w_kk_lb_n),
                      res_kk_lb$iterations, ok_kk_lb, note_kk_lb)
cat(sprintf("  %-22s wall=%7.4fs status=%d max_err=%.3e max_w=%.3f n_eff=%.1f\n",
            "leafblower_oris_soft", row_kk_lb$wall_s, res_kk_lb$status,
            row_kk_lb$max_error, row_kk_lb$max_w, row_kk_lb$n_eff))

# D-04 fallback: raking-with-acceleration is the documented path for
# zero-compression-benefit inputs -- this arm is the measurement backing
# that statement. Standing rule: raking's own canonical convergence
# (metric="kl", rule="improvement") is used untouched, no override; only
# `accelerate` and the shared iteration budget are set. Label by what the
# code does (accelerate=TRUE), not by an algorithm name ("SQUAREM" /
# "SRAA-m") the code does not itself assert -- resolving which acceleration
# scheme the flag selects is out of this phase's scope.
cat("  leafblower_raking_accelerated solving...\n")
kk_raking_call <- function() {
  harvest(df_kk, tgt_kk, method = "raking", accelerate = TRUE, bounds_mode = "unit",
          max_weight = 3, attach_weights = FALSE, max_iterations = kk_iter_budget)
}
bm_kk_raking  <- bench::mark(run = kk_raking_call(), iterations = 2, check = FALSE,
                              memory = FALSE, filter_gc = FALSE)
w_kk_raking   <- kk_raking_call()
res_kk_raking <- attr(w_kk_raking, "result")
w_kk_raking_n <- as.numeric(w_kk_raking)
max_error_kk_raking <- margin_max_error(w_kk_raking_n, df_kk, tgt_kk)
ok_kk_raking <- isTRUE(res_kk_raking$status %in% c(0L, 5L))
note_kk_raking <- sprintf(
  "accelerate=TRUE (D-04 raking fallback path for zero-compression-benefit inputs); raking's own canonical metric='kl',rule='improvement' default, no convergence override (standing rule); max_iterations=%d budget; status=%d, iterations=%d",
  kk_iter_budget, res_kk_raking$status, res_kk_raking$iterations)
row_kk_raking <- arm_row("known_limit_k20_uniform", n_kk, K_kk, cats_kk, m_cell_kk, 3,
                          "leafblower_raking_accelerated", as.numeric(bm_kk_raking$median),
                          max_error_kk_raking, max(w_kk_raking_n), min(w_kk_raking_n),
                          leafblower::design_effect(w_kk_raking_n),
                          leafblower::effective_sample_size(w_kk_raking_n),
                          res_kk_raking$iterations, ok_kk_raking, note_kk_raking)
cat(sprintf("  %-22s wall=%7.4fs status=%d max_err=%.3e max_w=%.3f n_eff=%.1f\n",
            "leafblower_raking_accelerated", row_kk_raking$wall_s, res_kk_raking$status,
            row_kk_raking$max_error, row_kk_raking$max_w, row_kk_raking$n_eff))

results <- rbind(results, row_kk_lb, row_kk_raking)

dir.create("benchmarks/results", showWarnings = FALSE, recursive = TRUE)
write.csv(results, "benchmarks/results/oris_soft_vs_competitors.csv", row.names = FALSE)
cat("\nWrote benchmarks/results/oris_soft_vs_competitors.csv\n")

# --- Machine and provenance capture (SC1: the figure must name its machine) ---
cpu_model <- tryCatch({
  if (file.exists("/proc/cpuinfo")) {
    lines <- readLines("/proc/cpuinfo")
    m <- grep("^model name", lines, value = TRUE)[1]
    if (!is.na(m)) trimws(sub("^model name\\s*:\\s*", "", m)) else NA_character_
  } else {
    NA_character_  # non-Linux host: degrade to NA rather than erroring
  }
}, error = function(e) NA_character_)

si <- sessionInfo()
# T-03-01: icarus and ReGenesees are deliberately absent from DESCRIPTION
# (D-09) and resolve from the ambient R library with no manifest pinning
# them -- this provenance line is the SOLE record of which build produced
# their published figures, so a substituted local package is visible as a
# provenance diff rather than a silent number change.
competitor_pkgs <- c("survey", "icarus", "ReGenesees", "optweight")
env_lines <- c(
  sprintf("R version: %s", R.version.string),
  sprintf("Platform: %s", R.version$platform),
  sprintf("BLAS: %s", si$BLAS),
  sprintf("LAPACK: %s", si$LAPACK),
  sprintf("CPU model: %s", if (is.na(cpu_model)) "NA (non-Linux host)" else cpu_model),
  sprintf("OMP_NUM_THREADS: %s", Sys.getenv("OMP_NUM_THREADS")),
  sprintf("OPENBLAS_NUM_THREADS: %s", Sys.getenv("OPENBLAS_NUM_THREADS")),
  sprintf("MKL_NUM_THREADS: %s", Sys.getenv("MKL_NUM_THREADS")),
  sprintf("leafblower: %s", as.character(utils::packageVersion("leafblower"))),
  vapply(competitor_pkgs, function(p)
    sprintf("%s: %s", p,
            if (requireNamespace(p, quietly = TRUE)) as.character(utils::packageVersion(p)) else "not installed"),
    character(1))
)
writeLines(env_lines, "benchmarks/results/oris_soft_vs_competitors_env.txt")
cat("Wrote benchmarks/results/oris_soft_vs_competitors_env.txt\n")
