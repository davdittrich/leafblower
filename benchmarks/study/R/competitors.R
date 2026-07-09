# benchmarks/study/R/competitors.R
#
# Uniform-contract adapters for 17 R competitor calibration/raking/balancing
# packages (registry.json R-arm ids). Ticket leafblower-2ouc.6 (WU-5).
#
# STRICT SEPARATION (user constraint, 2026-07-08): this file wraps COMPETITOR
# packages only. It never calls leafblower's C++ core, R/, or python/leafblower
# code -- leafblower's own adapters are WU-7, out of scope here.
#
# Every run_<registry_id>(problem) function returns the frozen contract v2
# shape (spec/contract.md Section 1):
#   list(weights_ref=, iterations=, status=, converged=, error_message=,
#        wall_time_s=, peak_rss_bytes=)
# `converged` is ALWAYS harness-recomputed here via common/metrics.R's
# margin_stats()$margin_linf <= problem$tol -- never a solver's self-report
# (contract.md Section 2.4). `status` is the adapter's own mapping of the
# wrapped package's native exit/exception into the 8-value harmonized enum
# (spec/status_enum.json); a package's own "success" can still coexist with
# converged==FALSE and vice versa (contract.md Section 2.3).
#
# Hyperparameters are sourced VERBATIM from spec/hyperparams.json (frozen,
# git-tagged before WU-REH rehearsal) -- they are not tuned per-problem, and a
# package's own defaults producing a poor margin match on a given problem is a
# reportable finding, not an adapter bug to paper over (see ebal below).

# ---- dependencies -----------------------------------------------------------
# Relative-to-repo-root paths (matches sibling benchmarks/study/R/
# leafblower_adapter.R's here() convention and CLAUDE.md's "run from repo
# root" invocation contract) -- deliberately NOT commandArgs()-derived, since
# that breaks under source()-from-a-different-script or testthat::test_file()
# (which chdir's into the test file's own directory).

# Guarded: run_worker() (run_arm.R:109-116) sources problem_io.R, then
# install_gen_resolver() monkey-patches its script-level .pio_resolve_data_ref
# so gen:instance_family?... data_refs resolve. An unconditional re-source
# here would redefine every problem_io.R function (including the patched
# .pio_resolve_data_ref) back to the stock, gen:-unaware version, silently
# breaking instance-family competitor cells (leafblower-2ouc.24). Skipping the
# source when load_problem_spec already exists preserves whatever loader
# (patched or stock) the caller already installed; standalone/test callers
# that haven't sourced problem_io.R yet still get it here.
if (!exists("load_problem_spec", mode = "function")) {
  source(file.path("benchmarks", "study", "common", "problem_io.R"))
}
if (!exists("margin_stats", mode = "function")) {
  source(file.path("benchmarks", "study", "common", "metrics.R"))
}

# ---- bad-arg signalling ------------------------------------------------------
# A distinct condition class so the harness wrapper (.comp_run) can map
# pre-solve input rejection to status="bad_arg" specifically, rather than the
# generic "error" bucket (status_enum.json: bad_arg = "Invalid argument /
# malformed problem input rejected before any solve iterate was produced").

.comp_bad_arg <- function(msg) {
  cond <- structure(
    class = c("comp_bad_arg_error", "error", "condition"),
    list(message = msg, call = sys.call(-1))
  )
  stop(cond)
}

# ---- peak RSS -----------------------------------------------------------

#' High-water-mark RSS in bytes (contract.md Section 2.7), Linux
#' /proc/self/status:VmHWM. Returns NA_real_ if unavailable (non-Linux).
.comp_peak_rss_bytes <- function() {
  path <- "/proc/self/status"
  if (!file.exists(path)) return(NA_real_)
  lines <- readLines(path, warn = FALSE)
  hwm <- grep("^VmHWM:", lines, value = TRUE)
  if (length(hwm) == 0) return(NA_real_)
  kb <- as.numeric(regmatches(hwm, regexpr("[0-9]+", hwm)))
  kb * 1024
}

# ---- weights_ref parquet convention ------------------------------------

#' weights/<solver>__<problem>__t<thread>__<build>.parquet (matches
#' common/ref_convex.R's ref_weights_path() convention; competitor rows are
#' always build="na" per contract.md Section 3 / registry.json).
.comp_weights_path <- function(solver_id, problem_id, thread = 1L, build = "na",
                                out_dir = "benchmarks/study/results") {
  file.path(out_dir, "weights",
            sprintf("%s__%s__t%d__%s.parquet", solver_id, problem_id, thread, build))
}

.comp_write_weights <- function(w, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  arrow::write_parquet(data.frame(weight = as.numeric(w)), path)
  invisible(path)
}

# ---- margin encodings ---------------------------------------------------
# Two conventions (see WU-5 investigation): FULL-dummy (all K categories per
# margin, no reference dropped) for packages whose solver tolerates/requires
# the resulting rank-deficient-but-consistent system (generalized inverse /
# IPF-style: sampling::calib, laeken::calibWeights, icarus::calibration,
# optweight, sbw); INTERCEPT + (K-1)-dummy classic regression encoding for
# packages that do a strict rcond() invertibility check (GECal, jointCalib's
# calib_el) or that natively build model.matrix(~margins) themselves
# (survey::calibrate, ReGenesees::e.calibrate) -- and for ebal, whose internal
# collinearity check rejects the full-dummy form even without an explicit
# intercept column (K-1 dummy, no intercept, is what ebal accepts).

.comp_groups <- function(problem) {
  setNames(lapply(problem$margins, function(m) as.character(problem$data[[m]])),
           problem$margins)
}

#' Full-dummy (no intercept) design matrix + concatenated population-count
#' targets, column-aligned. Used by sampling/laeken/icarus and as the base for
#' the ATT phantom-row trick (optweight/sbw/weightit-optweight).
.comp_fulldummy <- function(problem) {
  N <- sum(problem$design_weights)
  cols <- list()
  tot <- c()
  for (m in problem$margins) {
    f <- problem$data[[m]]
    Xi <- stats::model.matrix(~ f - 1)
    colnames(Xi) <- paste0(m, ".", levels(f))
    cols[[m]] <- Xi
    tgt <- problem$targets[[m]]
    names(tgt) <- paste0(m, ".", names(tgt))
    tot <- c(tot, tgt[colnames(Xi)] * N)
  }
  list(X = do.call(cbind, cols), totals = tot)
}

#' Intercept + (K-1)-dummy classic regression encoding + aligned totals
#' (intercept target = N; each non-reference dummy's target = its
#' population count). Used by GECal, jointCalib::calib_el.
.comp_regression <- function(problem) {
  N <- sum(problem$design_weights)
  df <- problem$data[, problem$margins, drop = FALSE]
  form <- stats::as.formula(paste("~", paste(problem$margins, collapse = " + ")))
  X <- stats::model.matrix(form, data = df)
  totals <- c(`(Intercept)` = N)
  for (m in problem$margins) {
    f <- problem$data[[m]]
    ref <- levels(f)[1]
    nonref <- setdiff(levels(f), ref)
    tgt <- problem$targets[[m]]
    nm <- paste0(m, nonref)
    totals <- c(totals, setNames(tgt[nonref] * N, nm))
  }
  totals <- totals[colnames(X)]
  list(X = X, totals = totals)
}

#' ATT "phantom target row" trick (full-dummy variant): append one synthetic
#' row per problem whose covariate values ARE the target proportions,
#' distinguished from the n real rows by a binary treat/focal indicator.
#' Turns single-sample calibration-to-population-totals into the ATT
#' balancing problem these packages are built for (ebal/optweight/sbw/
#' WeightIt were designed for causal-contrast balancing, not raw totals
#' calibration). Full-dummy variant: for optweight/sbw/WeightIt-optweight,
#' whose QP/conic solvers tolerate the resulting collinear-but-consistent
#' constraint system directly.
.comp_phantom <- function(problem) {
  fd <- .comp_fulldummy(problem)
  n <- nrow(problem$data)
  ph <- numeric(ncol(fd$X))
  names(ph) <- colnames(fd$X)
  for (m in problem$margins) {
    tgt <- problem$targets[[m]]
    nm <- paste0(m, ".", names(tgt))
    ph[nm] <- unname(tgt)
  }
  X_all <- rbind(fd$X, matrix(ph, nrow = 1, dimnames = list(NULL, names(ph))))
  list(X = fd$X, X_all = X_all, treat = c(rep(0, n), 1))
}

#' ATT phantom-row trick, K-1-dummy-no-intercept variant: ebal's internal
#' collinearity check rejects the full-dummy form (all K indicator columns
#' sum to a constant row, which ebal detects and refuses regardless of
#' whether an explicit intercept column is present); dropping one reference
#' level per margin avoids it. Used by ebal and WeightIt's ebal delegate.
.comp_phantom_kminus1 <- function(problem) {
  n <- nrow(problem$data)
  cols <- list()
  ph <- c()
  for (m in problem$margins) {
    f <- problem$data[[m]]
    lv <- levels(f)
    nonref <- lv[-1]
    Xi <- stats::model.matrix(~f)[, -1, drop = FALSE]
    colnames(Xi) <- paste0(m, ".", nonref)
    cols[[m]] <- Xi
    tgt <- problem$targets[[m]]
    ph <- c(ph, setNames(unname(tgt[nonref]), colnames(Xi)))
  }
  X <- do.call(cbind, cols)
  X_all <- rbind(X, matrix(ph, nrow = 1, dimnames = list(NULL, names(ph))))
  list(X = X, X_all = X_all, treat = c(rep(0, n), 1))
}

# ---- family -> native-method dispatch ------------------------------------
# icarus/laeken/ReGenesees/GECal(entropy)/WeightIt(delegate method) each wrap
# THREE (or two) distance families behind ONE registry id + a native
# method/entropy switch; the problem's own objective_families selects it.

.comp_family1 <- function(problem, adapter) {
  f <- problem$objective_families
  if (length(f) != 1) {
    .comp_bad_arg(sprintf(
      "%s: adapter requires exactly one objective_families entry, got %d (%s)",
      adapter, length(f), paste(f, collapse = ",")))
  }
  f[[1]]
}

.COMP_FAMILY_TO_RAKE_METHOD <- c(kl = "raking", chi2 = "linear", logit = "logit")

# ---- run() harness wrapper -----------------------------------------------
# Shared by every run_<id>(): timing, RSS, weights-parquet write (incl. the
# all-NaN hard-failure sentinel, contract.md Section 2.1), harness-recomputed
# `converged` (never a solver's self-report, contract.md Section 2.4), and the
# harness-side bound_violation upgrade for packages with no native bound
# check (status_enum.json: only applied when the package itself reported
# "converged" -- a native no_conv/error/infeasible classification is left as
# the more specific outcome).
#
# `solve_fn(problem)` must either throw (caught generically as "error", or as
# "bad_arg" via .comp_bad_arg()) or return
#   list(weights=<numeric(n)>, iterations=<int|NA>, status=<enum>,
#        error_message=<str|NULL>)

.comp_run <- function(solver_id, problem, solve_fn, out_dir = "benchmarks/study/results",
                       thread = 1L, build = "na") {
  n <- nrow(problem$data)
  t0 <- Sys.time()

  result <- tryCatch({
    r <- solve_fn(problem)
    if (length(r$weights) != n) {
      stop(sprintf("adapter '%s' returned weights of length %d, expected n=%d",
                    solver_id, length(r$weights), n))
    }
    r
  }, comp_bad_arg_error = function(e) {
    list(weights = rep(NaN, n), iterations = NA_integer_, status = "bad_arg",
         error_message = conditionMessage(e))
  }, error = function(e) {
    list(weights = rep(NaN, n), iterations = NA_integer_, status = "error",
         error_message = conditionMessage(e))
  })
  # contract.md Section 2.6: timer stops at weights-out (solve complete). The
  # margin_stats recompute + bound-violation upgrade + parquet write below are
  # post-run harness diagnostics (Section 2.4) and MUST NOT be inside the timer.
  wall_time_s <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  w <- as.numeric(result$weights)
  iterations <- result$iterations
  if (is.null(iterations) || length(iterations) == 0 || is.na(iterations)) {
    iterations <- NA_integer_
  } else {
    iterations <- as.integer(iterations)
  }
  status <- result$status
  error_message <- result$error_message

  hard_failure <- all(is.na(w))
  if (!hard_failure) {
    ms <- margin_stats(w, .comp_groups(problem), problem$targets)
    converged <- isTRUE(ms$margin_linf <= problem$tol)
    if (identical(status, "converged")) {
      eps <- 1e-8
      viol <- any(w < problem$bounds$min - eps, na.rm = TRUE) ||
        any(w > problem$bounds$max + eps, na.rm = TRUE)
      if (viol) status <- "bound_violation"
    }
  } else {
    converged <- FALSE
  }

  peak_rss_bytes <- .comp_peak_rss_bytes()

  weights_ref <- .comp_weights_path(solver_id, problem$id, thread = thread, build = build,
                                     out_dir = out_dir)
  .comp_write_weights(w, weights_ref)

  list(
    weights_ref = weights_ref,
    iterations = iterations,
    status = status,
    converged = converged,
    error_message = error_message,
    wall_time_s = wall_time_s,
    peak_rss_bytes = as.numeric(peak_rss_bytes)
  )
}

# =============================================================================
# survey::calibrate -- raking / linear / logit
# hyperparams.json: maxit=50, epsilon=1e-7, force=FALSE; bounds required
# (and problem-supplied) only for calfun="logit".
# =============================================================================

.comp_survey_solve <- function(problem, calfun) {
  df <- problem$data
  df$.comp_w <- problem$design_weights
  des <- survey::svydesign(ids = ~1, weights = ~.comp_w, data = df)
  reg <- .comp_regression(problem)
  form <- stats::as.formula(paste("~", paste(problem$margins, collapse = " + ")))

  args <- list(design = des, formula = form, population = reg$totals,
               calfun = calfun, maxit = 50, epsilon = 1e-7, force = FALSE)
  if (calfun == "logit") {
    if (!is.finite(problem$bounds$min) || !is.finite(problem$bounds$max)) {
      .comp_bad_arg("survey::calibrate calfun='logit' requires finite bounds on both sides")
    }
    args$bounds <- c(problem$bounds$min, problem$bounds$max)
  }

  warned <- character(0)
  cal <- withCallingHandlers(
    do.call(survey::calibrate, args),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  w <- as.numeric(stats::weights(cal))
  no_conv <- any(grepl("did not converge|failed to converge", warned, ignore.case = TRUE))
  list(weights = w, iterations = NA_integer_,
       status = if (no_conv) "no_conv" else "converged",
       error_message = if (length(warned)) paste(warned, collapse = "; ") else NULL)
}

run_survey_calibrate_raking <- function(problem) {
  .comp_run("survey_calibrate_raking", problem, function(p) .comp_survey_solve(p, "raking"))
}
run_survey_calibrate_linear <- function(problem) {
  .comp_run("survey_calibrate_linear", problem, function(p) .comp_survey_solve(p, "linear"))
}
run_survey_calibrate_logit <- function(problem) {
  .comp_run("survey_calibrate_logit", problem, function(p) .comp_survey_solve(p, "logit"))
}

# =============================================================================
# sampling::calib -- linear / logit
# hyperparams.json: max_iter=500; bounds c(low,upp) required (problem-supplied)
# only for method="logit". NA in the returned g-weight vector signals
# non-convergence (package convention, not an exception).
# =============================================================================

.comp_sampling_solve <- function(problem, method) {
  fd <- .comp_fulldummy(problem)
  args <- list(Xs = fd$X, d = problem$design_weights, total = fd$totals,
               method = method, max_iter = 500)
  if (method == "logit") {
    if (!is.finite(problem$bounds$min) || !is.finite(problem$bounds$max)) {
      .comp_bad_arg("sampling::calib method='logit' requires finite bounds on both sides")
    }
    args$bounds <- c(low = problem$bounds$min, upp = problem$bounds$max)
  }
  g <- do.call(sampling::calib, args)
  if (any(is.na(g))) {
    return(list(weights = rep(NA_real_, length(g)), iterations = NA_integer_,
                status = "no_conv",
                error_message = "sampling::calib returned NA g-weight(s) (package's own non-convergence signal)"))
  }
  list(weights = problem$design_weights * g, iterations = NA_integer_,
       status = "converged", error_message = NULL)
}

run_sampling_calib_linear <- function(problem) {
  .comp_run("sampling_calib_linear", problem, function(p) .comp_sampling_solve(p, "linear"))
}
run_sampling_calib_logit <- function(problem) {
  .comp_run("sampling_calib_logit", problem, function(p) .comp_sampling_solve(p, "logit"))
}

# =============================================================================
# anesrake::anesrake
# hyperparams.json documented defaults verbatim: cap=5, maxit=1000,
# type="pctlim", pctlim=5, nlim=5, choosemethod="total", iterate=TRUE,
# convcrit=0.01, force1=TRUE. Home-turf golden: bundled anes04 dataset (a
# synthetic toy sample already within pctlim=5% of its own targets makes
# anesrake refuse to run at all -- "No variables are off by more than 5
# percent" -- which is the package's documented behaviour, not a bug).
# =============================================================================

.comp_anesrake_solve <- function(problem) {
  df <- as.data.frame(lapply(problem$data[problem$margins], as.character),
                       stringsAsFactors = FALSE)
  names(df) <- problem$margins
  # anesrake::anesrake()'s OWN source does `dataframe <- dataframe[filter ==
  # 1, ]` (default scalar filter=1) BEFORE its names(inputter) %in%
  # names(dataframe) guard. For a single-margin (K=1) data.frame this is R's
  # ordinary `[.data.frame` two-index drop=TRUE default silently collapsing
  # the 1-column result to a plain unnamed vector -- names(dataframe) then
  # returns NULL and anesrake's own guard fails with a misleading "variable
  # not found in the data frame" error, even though the raking variable IS
  # present. A harmless placeholder column (excluded from the target list,
  # so never raked on) keeps ncol(df) >= 2 and avoids the drop.
  df$.comp_placeholder <- 1
  target <- lapply(problem$margins, function(m) problem$targets[[m]])
  names(target) <- problem$margins
  caseid <- seq_len(nrow(df))

  res <- tryCatch({
    withCallingHandlers(
      anesrake::anesrake(target, df, caseid = caseid, weightvec = problem$design_weights,
                          cap = 5, verbose = FALSE, maxit = 1000, type = "pctlim",
                          pctlim = 5, nlim = 5, choosemethod = "total", iterate = TRUE,
                          convcrit = 0.01, force1 = TRUE),
      message = function(m) invokeRestart("muffleMessage")
    )
  }, error = function(e) {
    if (grepl("off by more than", conditionMessage(e), fixed = TRUE)) {
      structure(list(msg = conditionMessage(e)), class = "comp_anesrake_noop")
    } else {
      stop(e)
    }
  })

  if (inherits(res, "comp_anesrake_noop")) {
    # anesrake's own gate found nothing exceeding pctlim=5%: the sample
    # already matches its targets closely enough that the package refuses to
    # rake at all. Design weights unchanged is the correct, sum-preserving,
    # already-at-target result -- not a failure (contract.md Section 2.5: a
    # non-null error_message does not itself imply status != "converged").
    return(list(weights = problem$design_weights, iterations = NA_integer_,
                status = "converged", error_message = res$msg))
  }

  w <- as.numeric(res$weightvec)
  it <- suppressWarnings(as.integer(res$iterations))
  status <- if (grepl("Complete convergence", res$converge, fixed = TRUE)) "converged" else "no_conv"
  list(weights = w, iterations = if (length(it) && !is.na(it)) it else NA_integer_,
       status = status, error_message = NULL)
}

run_anesrake <- function(problem) .comp_run("anesrake", problem, .comp_anesrake_solve)

# =============================================================================
# ipfr::ipu
# hyperparams.json: max_iterations=30, min_weight=0.0001 (installed API's
# corresponding formal parameter is named weight_floor), secondary_importance
# not used (no secondary targets in this contract), min_ratio/max_ratio null
# (package defaults max_ratio=10000/min_ratio=0.0001 apply).
# =============================================================================

.comp_ipfr_solve <- function(problem) {
  N <- sum(problem$design_weights)
  seed <- problem$data[problem$margins]
  seed$id <- seq_len(nrow(seed))
  seed$weight <- problem$design_weights
  primary_targets <- lapply(problem$margins, function(m) {
    tgt <- as.list(problem$targets[[m]] * N)
    as.data.frame(tgt)
  })
  names(primary_targets) <- problem$margins

  warned <- character(0)
  res <- withCallingHandlers(
    ipfr::ipu(primary_seed = seed, primary_targets = primary_targets,
              max_iterations = 30, weight_floor = 1e-4, relative_gap = 0.01,
              absolute_diff = 10),
    warning = function(w) { warned <<- c(warned, conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  w <- as.numeric(res$weight_tbl$weight)
  no_conv <- any(grepl("did not converge|max.*iteration", warned, ignore.case = TRUE))
  list(weights = w, iterations = NA_integer_,
       status = if (no_conv) "no_conv" else "converged",
       error_message = if (length(warned)) paste(warned, collapse = "; ") else NULL)
}

run_ipfr <- function(problem) .comp_run("ipfr", problem, .comp_ipfr_solve)

# =============================================================================
# ebal::ebalance
# hyperparams.json documented defaults verbatim: max.iterations=200,
# constraint.tolerance=1 (RAW COUNT-scale absolute moment imbalance, not a
# proportion -- on small-n problems this frozen default yields a genuinely
# loose margin match; that is the package's documented behaviour and is
# reported via converged==FALSE by the harness recompute, not tightened here).
# ATT phantom-row trick, K-1-dummy-no-intercept encoding (ebal's collinearity
# check rejects the full-dummy form).
# =============================================================================

.comp_ebal_solve <- function(problem) {
  ph <- .comp_phantom_kminus1(problem)
  res <- ebal::ebalance(Treatment = ph$treat, X = ph$X_all,
                         max.iterations = 200, constraint.tolerance = 1)
  w <- as.numeric(res$w)
  status <- if (isTRUE(res$converged)) "converged" else "no_conv"
  list(weights = w, iterations = NA_integer_, status = status, error_message = NULL)
}

run_ebal <- function(problem) .comp_run("ebal", problem, .comp_ebal_solve)

# =============================================================================
# optweight::optweight -- entropy (kl family) / linf (minimax family)
# hyperparams.json: tols=0 (exact balance), min.w=1e-8. ATT phantom-row
# trick, full-dummy encoding (optweight's QP/conic solver tolerates the
# collinear-but-consistent system directly, unlike ebal).
# =============================================================================

.comp_optweight_solve <- function(problem, norm) {
  ph <- .comp_phantom(problem)
  treat <- factor(ifelse(ph$treat == 1, "target", "sample"), levels = c("sample", "target"))
  dfall <- data.frame(treat = treat, ph$X_all)
  form <- stats::as.formula(paste("treat ~", paste(colnames(ph$X_all), collapse = " + ")))
  fit <- optweight::optweight(form, data = dfall, tols = 0, estimand = "ATT",
                               focal = "target", norm = norm, min.w = 1e-8)
  w <- as.numeric(fit$weights[treat == "sample"])

  if (norm == "entropy") {
    st <- fit$info$status
    it <- suppressWarnings(as.integer(fit$info$iter))
    status <- if (identical(st, "solved")) "converged"
    else if (grepl("infeasible", st, fixed = TRUE)) "infeasible"
    else "no_conv"
    msg <- if (status != "converged") paste0("optweight (SCS) status: ", st) else NULL
  } else {
    st <- fit$info$status_message
    it <- suppressWarnings(as.integer(fit$info$info$simplex_iteration_count))
    status <- if (identical(st, "Optimal")) "converged"
    else if (grepl("infeasible", st, ignore.case = TRUE)) "infeasible"
    else "no_conv"
    msg <- if (status != "converged") paste0("optweight (HiGHS) status: ", st) else NULL
  }
  list(weights = w, iterations = if (length(it) && !is.na(it)) it else NA_integer_,
       status = status, error_message = msg)
}

run_optweight_entropy <- function(problem) {
  .comp_run("optweight_entropy", problem, function(p) .comp_optweight_solve(p, "entropy"))
}
run_optweight_linf <- function(problem) {
  .comp_run("optweight_linf", problem, function(p) .comp_optweight_solve(p, "linf"))
}

# =============================================================================
# WeightIt::weightit -- delegates to ebal (kl family) / optweight (minimax
# family) per hyperparams.json note ("Driver must record which method
# WeightIt is invoked with and inherit that engine's frozen entry").
# =============================================================================

.comp_weightit_solve <- function(problem) {
  fam <- .comp_family1(problem, "weightit")
  method <- unname(c(kl = "ebal", minimax = "optweight")[fam])
  if (is.na(method)) {
    .comp_bad_arg(sprintf("weightit: unsupported objective family '%s' (expected kl|minimax)", fam))
  }
  ph <- if (method == "ebal") .comp_phantom_kminus1(problem) else .comp_phantom(problem)
  dfall <- data.frame(treat = ph$treat, ph$X_all)
  form <- stats::as.formula(paste("treat ~", paste(colnames(ph$X_all), collapse = " + ")))
  fit <- WeightIt::weightit(form, data = dfall, method = method, estimand = "ATT")
  w <- as.numeric(fit$weights[ph$treat == 0])
  list(weights = w, iterations = NA_integer_, status = "converged", error_message = NULL)
}

run_weightit <- function(problem) .comp_run("weightit", problem, .comp_weightit_solve)

# =============================================================================
# ReGenesees::e.calibrate -- linear / raking / logit (family dispatch, as
# icarus/laeken below). hyperparams.json: maxit=50, epsilon=1e-7, force=TRUE.
# =============================================================================

.comp_regenesees_solve <- function(problem, calfun) {
  df <- problem$data
  df$.comp_id <- seq_len(nrow(df))
  df$.comp_w <- problem$design_weights
  des <- ReGenesees::e.svydesign(data = df, ids = ~.comp_id, weights = ~.comp_w,
                                  check.data = TRUE)
  form <- stats::as.formula(paste("~", paste(problem$margins, collapse = " + "), "- 1"))
  poptot <- ReGenesees::pop.template(des, calmodel = form, partition = FALSE)
  N <- sum(problem$design_weights)
  for (cn in colnames(poptot)) {
    m <- problem$margins[vapply(problem$margins, function(mm) startsWith(cn, mm), logical(1))][1]
    lvl <- substring(cn, nchar(m) + 1)
    poptot[1, cn] <- unname(problem$targets[[m]][lvl]) * N
  }
  args <- list(design = des, df.population = poptot, calmodel = form, partition = FALSE,
               calfun = calfun, maxit = 50, epsilon = 1e-7, force = TRUE)
  if (calfun == "logit") {
    if (!is.finite(problem$bounds$min) || !is.finite(problem$bounds$max)) {
      .comp_bad_arg("ReGenesees::e.calibrate calfun='logit' requires finite bounds on both sides")
    }
    args$bounds <- c(problem$bounds$min, problem$bounds$max)
  }
  cal <- do.call(ReGenesees::e.calibrate, args)
  w <- as.numeric(stats::weights(cal))
  list(weights = w, iterations = NA_integer_, status = "converged", error_message = NULL)
}

run_regenesees <- function(problem) {
  .comp_run("regenesees", problem, function(p) {
    fam <- .comp_family1(p, "regenesees")
    calfun <- unname(.COMP_FAMILY_TO_RAKE_METHOD[fam])
    if (is.na(calfun)) .comp_bad_arg(sprintf("regenesees: unsupported objective family '%s'", fam))
    .comp_regenesees_solve(p, calfun)
  })
}

# =============================================================================
# icarus::calibration -- linear / raking / logit (family dispatch).
# hyperparams.json: maxIter=2500, calibTolerance=1e-6; bounds required
# (problem-supplied) only for method="logit".
# =============================================================================

.comp_icarus_solve <- function(problem, method) {
  N <- sum(problem$design_weights)
  df <- problem$data[, problem$margins, drop = FALSE]
  df$.comp_w <- problem$design_weights
  mm <- icarus::newMarginMatrix()
  for (m in problem$margins) {
    lv <- levels(problem$data[[m]])
    mm <- icarus::addMargin(mm, m, as.numeric(problem$targets[[m]][lv]))
  }
  args <- list(data = df, marginMatrix = mm, colWeights = ".comp_w", method = method,
               popTotal = N, pct = TRUE, maxIter = 2500, calibTolerance = 1e-6,
               description = FALSE)
  if (method == "logit") {
    if (!is.finite(problem$bounds$min) || !is.finite(problem$bounds$max)) {
      .comp_bad_arg("icarus::calibration method='logit' requires finite bounds on both sides")
    }
    args$bounds <- c(low = problem$bounds$min, upp = problem$bounds$max)
  }
  w <- as.numeric(do.call(icarus::calibration, args))
  list(weights = w, iterations = NA_integer_, status = "converged", error_message = NULL)
}

run_icarus <- function(problem) {
  .comp_run("icarus", problem, function(p) {
    fam <- .comp_family1(p, "icarus")
    method <- unname(.COMP_FAMILY_TO_RAKE_METHOD[fam])
    if (is.na(method)) .comp_bad_arg(sprintf("icarus: unsupported objective family '%s'", fam))
    .comp_icarus_solve(p, method)
  })
}

# =============================================================================
# laeken::calibWeights -- raking / linear / logit (family dispatch).
# hyperparams.json: maxit=500, tol=1e-6, eps=.Machine$double.eps,
# bounds=c(0,10) package default (overridden with problem bounds for logit).
# Full-dummy encoding (tolerated via generalized inverse, like sampling::calib).
# =============================================================================

.comp_laeken_solve <- function(problem, method) {
  fd <- .comp_fulldummy(problem)
  args <- list(X = fd$X, d = problem$design_weights, totals = fd$totals, method = method,
               maxit = 500, tol = 1e-6, eps = .Machine$double.eps)
  if (method == "logit") {
    if (!is.finite(problem$bounds$min) || !is.finite(problem$bounds$max)) {
      .comp_bad_arg("laeken::calibWeights method='logit' requires finite bounds on both sides")
    }
    args$bounds <- c(problem$bounds$min, problem$bounds$max)
  } else {
    args$bounds <- c(0, 10)
  }
  g <- do.call(laeken::calibWeights, args)
  if (any(is.na(g))) {
    return(list(weights = rep(NA_real_, length(g)), iterations = NA_integer_,
                status = "no_conv",
                error_message = "laeken::calibWeights returned NA g-weight(s) (package's own non-convergence signal)"))
  }
  list(weights = problem$design_weights * g, iterations = NA_integer_,
       status = "converged", error_message = NULL)
}

run_laeken <- function(problem) {
  .comp_run("laeken", problem, function(p) {
    fam <- .comp_family1(p, "laeken")
    method <- unname(.COMP_FAMILY_TO_RAKE_METHOD[fam])
    if (is.na(method)) .comp_bad_arg(sprintf("laeken: unsupported objective family '%s'", fam))
    .comp_laeken_solve(p, method)
  })
}

# =============================================================================
# GECal::GEcalib -- family -> entropy dispatch (kl->ET exp-tilting/raking,
# chi2->SL squared-loss, logit->CE), method GEC0 (no debiasing covariate) when
# design weights are constant, GEC (+ g(dweight) debiasing term, population
# total unknown -> NA per package's own documented example) when not.
# hyperparams.json: xtol=1e-16, maxit=1e5, weight.scale=1, G.scale=1,
# allowSingular=TRUE. Intercept+(K-1)-dummy regression encoding (GEcalib does
# a strict rcond() invertibility check that rejects the full-dummy form).
# =============================================================================

.comp_gecal_solve <- function(problem) {
  fam <- .comp_family1(problem, "gecal")
  entropy <- unname(c(kl = "ET", chi2 = "SL", logit = "CE")[fam])
  if (is.na(entropy)) .comp_bad_arg(sprintf("gecal: unsupported objective family '%s'", fam))

  reg <- .comp_regression(problem)
  data <- problem$data[, problem$margins, drop = FALSE]
  .d <- problem$design_weights
  main_form <- paste("~", paste(problem$margins, collapse = " + "))

  if (length(unique(.d)) <= 1L) {
    fit <- GECal::GEcalib(stats::as.formula(main_form), dweight = .d, data = data,
                           const = reg$totals, method = "GEC0", entropy = entropy,
                           xtol = 1e-16, maxit = 100000, weight.scale = 1, G.scale = 1,
                           allowSingular = TRUE)
  } else {
    form <- stats::as.formula(paste(main_form, "+ g(.d)"))
    const <- c(reg$totals, NA_real_)
    names(const)[length(const)] <- "g(.d)"
    fit <- GECal::GEcalib(form, dweight = .d, data = data, const = const,
                           method = "GEC", entropy = entropy, xtol = 1e-16, maxit = 100000,
                           weight.scale = 1, G.scale = 1, allowSingular = TRUE)
  }
  w <- as.numeric(fit$w)
  list(weights = w, iterations = NA_integer_, status = "converged", error_message = NULL)
}

run_gecal <- function(problem) .comp_run("gecal", problem, .comp_gecal_solve)

# =============================================================================
# jointCalib -- calib_el() directly, NOT joint_calib().
#
# BLOCKER RESOLVED: joint_calib() unconditionally rejects totals-only
# calibration --
#   stop("The `formula_quantiles` parameter is required. If you want to use
#         standard calibration, we recommend using the `survey`, `sampling`,
#         `laeken` or `ebal` packages.")
# (verified via jointCalib:::joint_calib source, no argument combination
# bypasses this gate). jointCalib exports its lower-level empirical-likelihood
# calibration ENGINE directly as calib_el(X, d, totals, maxit, tol, eps) --
# this is exactly the totals-only calibration primitive joint_calib() itself
# calls internally, with no formula_quantiles requirement. Uses the same
# intercept+(K-1)-dummy encoding as GECal (calib_el does X[,-1] internally,
# i.e. it expects an explicit intercept as column 1).
# hyperparams.json control_calib defaults: maxit=50, tol=1e-8,
# eps=.Machine$double.eps.
# =============================================================================

.comp_jointcalib_solve <- function(problem) {
  reg <- .comp_regression(problem)
  msgs <- character(0)
  g <- withCallingHandlers(
    jointCalib::calib_el(X = reg$X, d = problem$design_weights, totals = reg$totals,
                          maxit = 50, tol = 1e-8, eps = .Machine$double.eps),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )
  w <- problem$design_weights * g
  status <- if (length(msgs)) "no_conv" else "converged"
  list(weights = w, iterations = NA_integer_, status = status,
       error_message = if (length(msgs)) paste(msgs, collapse = "; ") else NULL)
}

run_jointcalib <- function(problem) .comp_run("jointcalib", problem, .comp_jointcalib_solve)

# =============================================================================
# sbw::sbw
# hyperparams.json documented defaults verbatim: bal_alg=TRUE (algorithmic
# grid search over bal_gri), bal_std="group", sol_nam="quadprog",
# wei_sum=TRUE, wei_pos=TRUE, par_est="att" -- i.e. the SAME ATT phantom-row
# trick as ebal/optweight (par_est="aux", a direct totals-only path, exists
# but is NOT the frozen hyperparams choice). Full-dummy encoding (sbw's QP
# tolerates it, like optweight). sol_nam="quadprog" reports status=NA
# unconditionally (package's own documented limitation: "For solver
# 'quadprog', status code is missing, therefore status=NA").
# =============================================================================

.comp_sbw_solve <- function(problem) {
  ph <- .comp_phantom(problem)
  dat <- data.frame(treat = ph$treat, ph$X_all)
  out <- sbw::sbw(dat = dat, ind = "treat",
                   bal = list(bal_cov = colnames(ph$X_all), bal_alg = TRUE, bal_std = "group"),
                   sol = list(sol_nam = "quadprog", sol_dis = FALSE),
                   par = list(par_est = "att"), mes = FALSE)
  w <- as.numeric(out$dat_weights$sbw_weights[ph$treat == 0])
  list(weights = w, iterations = NA_integer_, status = "converged", error_message = NULL)
}

run_sbw <- function(problem) .comp_run("sbw", problem, .comp_sbw_solve)

# =============================================================================
# nonprobsvy::nonprob -- GEE calibration-constraint dual (Chen, Li & Wu 2021
# corrected-estimating-equations calibration; est_method="gee" solved via
# nleqslv with its PACKAGE-DEFAULT global strategy "dbldog", double-dogleg
# trust region). This is DESIGN.md Section 2 line 71's "double-dogleg TR
# dual" -- the closest analog to leafblower's newton_kl solver.
#
# hyperparams.json's frozen nonprobsvy entry records nonprob()'s raw
# control_sel() SIGNATURE DEFAULT (est_method="mle", plain Newton-Raphson MLE
# on the corrected log-likelihood) in its `params` block, but its own
# convergence_bounds_quantity field identifies est_method="gee" as "the
# Chen-Li-Wu calibration-constraint dual (closest Newton-KL analog)" --
# matching this competitor's registry slot exactly, and est_method="mle"
# would not be TR-dual at all. Following the same driver-selects-the-analog
# precedent already used for the WeightIt entry ("Driver must record which
# method WeightIt is invoked with and inherit that engine's frozen entry"),
# this adapter invokes est_method="gee" and leaves every other control_sel()
# knob at ITS OWN package default -- nleqslv_method="Broyden",
# nleqslv_global="dbldog", gee_h_fun=1, epsilon=1e-4, maxit=500, all verbatim
# from hyperparams.json's recorded values (zero per-problem tuning).
#
# API-fit finding (WU-5b investigation): nonprob() is architected for
# non-probability-sample POPULATION-MEAN estimation and unconditionally
# requires a `target` formula (the response variable of interest) -- even a
# pure calibration-weight extraction with no outcome model errors with
# "Please provide the `target` argument as a formula." if `target` is
# omitted. Verified empirically that this is a mandatory-argument formality
# ONLY: the returned per-observation weights (`weights.nonprob()` S3 method)
# are bit-identical across two independent, unrelated placeholder `target`
# columns -- the GEE/MLE calibration-constraint weight construction depends
# only on the `selection` covariates and `pop_totals`, never on the target
# response. A placeholder numeric column therefore satisfies the API gate
# without influencing the returned weights: a genuine calibration-weight-out
# path exists, not a forced mapping. Intercept+(K-1)-dummy regression
# encoding (matches GECal/jointCalib: nonprob()'s internal
# model.matrix(selection, data) needs pop_totals named to its own colnames).
#
# No native bound-constraint argument exists on nonprob() (confirmed via
# args(nonprob) -- no bounds= parameter): relies on .comp_run()'s
# harness-side bound_violation upgrade, like GECal/jointCalib/WeightIt/sbw.
#
# Numerical-robustness finding: on a large, poorly-scaled real dataset
# (nonprobsvy's own bundled admin/jvs job-vacancy pair, n=9344) this same
# invocation produces nleqslv's own "Ill-conditioned Jacobian"/"Jacobian is
# singular" warnings and returns weights that do NOT hit the declared
# pop_totals -- caught correctly by the harness-recomputed `converged`
# (never trusting the package's own non-error exit), not silently reported
# as success. The home-turf golden below uses the well-conditioned K=3
# apistrat fixture (shared with the survey::calibrate goldens above), where
# the GEE solve is exact to machine precision.
# =============================================================================

.comp_nonprobsvy_solve <- function(problem) {
  reg <- .comp_regression(problem)
  df <- problem$data[, problem$margins, drop = FALSE]
  # Mandatory API formality only -- verified not to affect the returned
  # weights (see comment block above).
  df$.comp_target <- seq_len(nrow(df))
  sel_form <- stats::as.formula(paste("~", paste(problem$margins, collapse = " + ")))

  fit <- nonprobsvy::nonprob(
    data = df, selection = sel_form, target = ~ .comp_target,
    pop_totals = reg$totals, method_selection = "logit",
    control_selection = nonprobsvy::control_sel(
      est_method = "gee", gee_h_fun = 1, epsilon = 1e-4, maxit = 500),
    se = FALSE, verbose = FALSE)
  w <- as.numeric(stats::weights(fit))
  list(weights = w, iterations = NA_integer_, status = "converged", error_message = NULL)
}

run_nonprobsvy <- function(problem) .comp_run("nonprobsvy", problem, .comp_nonprobsvy_solve)
