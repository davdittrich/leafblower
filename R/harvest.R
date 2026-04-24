#' Generate calibrated weights (drop-in for autumn::harvest)
#'
#' @param data A data frame containing columns to calibrate.
#' @param target A named list of named numeric vectors (variable -> proportions).
#' @param min_weight Lower bound on weights. Default 0 (no lower bound).
#' @param max_weight Upper bound on weights. Default 5.
#' @param method One of "auto", "ieppa", "lbfgsb", "raking".
#'   \itemize{
#'     \item \code{"ieppa"}: paper-faithful algBCD at C=0 (cell-compressed
#'       Sinkhorn, Chu-Liang-Toh-Yang 2022 arXiv:2011.14312).
#'     \item \code{"raking"}: classical IPF + Dykstra box + hyperplane
#'       projections (Deming-Stephan 1940 / Csiszar 1975 + Boyle-Dykstra
#'       1986). Renamed from the prior misnamed "iEPPA" hybrid.
#'     \item \code{"lbfgsb"}: L-BFGS-B on the Deville-Sarndal logit dual.
#'     \item \code{"auto"} (default): currently routes unconditionally to
#'       \code{"ieppa"}. Benchmark-driven routing refinement is planned.
#'   }
#' @param verbose Integer verbosity: 0=silent, 1=progress, 2=debug.
#' @param max_iterations Maximum inner BCD iterations per outer step. Default 500.
#' @param start_weights Starting weights vector or NULL (uniform).
#' @param attach_weights If TRUE, return data frame with weights column. Default TRUE.
#' @param weight_column Name of weight column. Default "weights".
#' @param convergence Named list/vector; "absolute" maps to tol_abs.
#' @param bounds_mode One of "cell" (default) or "unit". Controls per-observation
#'   vs cell-aggregate bound enforcement.
#' @param homotopy_levels Number of homotopy levels (default 1 = disabled). Values > 1
#'   progressively tighten max_weight from homotopy_start_factor to homotopy_end_factor
#'   across n levels.
#' @param homotopy_start_factor Starting max_weight multiplier (default 1.0 = no change).
#' @param homotopy_end_factor Ending max_weight multiplier (default 1.0 = no change).
#' @param homotopy_budget_p Budget split fraction across homotopy levels (default 0.5).
#' @param scheduler Margin sweep scheduler: "round_robin" (default) or "greedy"
#'   (Greenkhorn priority).
#' @param eta_schedule ALM penalty schedule: "fixed" (default) or "tang_dynamic".
#' @param eta_start Starting ALM penalty multiplier (default 1.0).
#' @param eta_end Ending ALM penalty multiplier (default 1.0).
#' @param eta_schedule_power Power for Tang-eta schedule interpolation (default 0.5).
#' @param select_params Ignored with verbose >= 2 note.
#' @param select_function Ignored.
#' @param error_function Ignored.
#' @param adaptive_order Ignored.
#' @param enforce_mean Ignored (retained for compatibility).
#' @param accelerate Ignored.
#' @param add_na_proportion Not supported in v1; raises error if TRUE.
#' @param auto_collapse Not supported in v1; raises error if TRUE.
#' @param collapse_vars Not supported in v1; raises error if TRUE.
#' @param target_map Passed through for data-frame target format handling.
#' @param ... Additional arguments ignored.
#' @return data frame with weights column if attach_weights=TRUE, else numeric vector.
#' @export
harvest <- function(
  data,
  target,
  min_weight       = 0,
  max_weight       = 5,
  method           = "ieppa",
  verbose          = 0,
  max_iterations   = 500,
  start_weights    = NULL,
  attach_weights   = TRUE,
  weight_column    = "weights",
  convergence      = list(),
  bounds_mode      = "cell",
  # --- new overlay knobs (all default off / identity) ---
  homotopy_levels       = 1L,
  homotopy_start_factor = 1.0,
  homotopy_end_factor   = 1.0,
  homotopy_budget_p     = 0.5,
  scheduler             = c("round_robin", "greedy"),
  eta_schedule          = c("fixed", "tang_dynamic"),
  eta_start             = 1.0,
  eta_end               = 1.0,
  eta_schedule_power    = 0.5,
  # --- end new ---
  select_params    = NULL,
  select_function  = NULL,
  error_function   = NULL,
  adaptive_order   = NULL,
  enforce_mean     = TRUE,
  accelerate       = FALSE,
  add_na_proportion = FALSE,
  auto_collapse    = FALSE,
  collapse_vars    = NULL,
  target_map       = NULL,
  design_weights   = NULL,
  ...
) {
  # Not-in-v1 hard stops
  if (!isFALSE(add_na_proportion))
    stop("add_na_proportion is not supported in leafblower v1.")
  if (isTRUE(auto_collapse))
    stop("auto_collapse is not supported in leafblower v1.")
  if (!is.null(collapse_vars))
    stop("collapse_vars is not supported in leafblower v1.")

  target  <- parse_target(target, target_map)
  method  <- map_method(method, verbose)
  # Emit deprecation warning if convergence['pct'] supplied.
  # tol_abs is forwarded to C_rk_calibrate as the 9th argument.
  tol_abs <- parse_convergence(convergence)

  # design_weights: used as start_weights when supplied (normalized to mean=1 by normalize_start_weights)
  if (!is.null(design_weights) && is.null(start_weights)) {
    start_weights <- design_weights
  }
  sw_vec  <- normalize_start_weights(start_weights, nrow(data))

  bounds_mode_char <- parse_bounds_mode(bounds_mode)
  bounds_mode_int  <- match(bounds_mode_char, c("cell", "unit")) - 1L

  # Overlay arg resolution
  scheduler    <- match.arg(scheduler)
  eta_schedule <- match.arg(eta_schedule)

  # Ignored-param verbose notes
  # enforce_mean is always TRUE: normalization is unconditional (line ~86).
  ignored <- c("select_params", "select_function", "error_function",
                "adaptive_order", "accelerate", "enforce_mean")
  supplied_ignored <- intersect(ignored, names(match.call(expand.dots = FALSE)))
  if (verbose >= 2 && length(supplied_ignored) > 0)
    message("[leafblower] Ignoring autumn params: ", paste(supplied_ignored, collapse = ", "))

  raw <- .Call("C_rk_calibrate",
               data,
               target,
               as.double(min_weight),
               as.double(max_weight),
               as.character(method),
               as.integer(verbose),
               as.integer(max_iterations),
               sw_vec,
               as.double(tol_abs),
               as.integer(bounds_mode_int),
               as.integer(homotopy_levels),
               as.double(homotopy_start_factor),
               as.double(homotopy_end_factor),
               as.double(homotopy_budget_p),
               as.character(scheduler),
               as.character(eta_schedule),
               as.double(eta_start),
               as.double(eta_end),
               as.double(eta_schedule_power),
               PACKAGE = "leafblower")

  weights <- raw$weights
  calib_result    <- raw$result

  # Check hard-stop statuses before normalization: status 2/3 mean weights are
  # meaningless; normalizing near-zero weights before stopping produces NaN.
  if (calib_result$status == 2L)
    stop("leafblower: infeasible problem \u2014 persistent empty cell with positive target (detected after 5 consecutive outer iterations).")
  if (calib_result$status == 3L)
    stop("leafblower: invalid arguments \u2014 ", calib_result$message)

  # Solver returns sum(weights) = n (enforced in src/ieppa.cpp, src/raking.cpp,
  # src/lbfgsb_solver.cpp per user directive 2026-04-24). No wrapper-level
  # normalization — removing it preserves the bounds_mode="unit" strict-bounds
  # guarantee (ieppa's water-fill clamps are final; not re-pushed by post-scale).

  # NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping
  # here would break sum(weights * d) == target totals when per-cell mixing
  # parameters d are non-uniform: individual weights may legitimately exceed
  # per-cell bounds after expansion even when cell aggregates are in range.
  # The iEPPA/LBFGSB solvers enforce bounds on the cell aggregate X[c], which
  # is the invariant that preserves calibration. See
  # tests/testthat/test-ieppa-nonuniform-d.R.

  if (calib_result$status == 1L)
    warning("leafblower: calibration did not converge (max_error=",
            signif(calib_result$max_error, 3), "). Weights reflect last iterate.")

  if (!is.null(calib_result$n_bounds_violated) && calib_result$n_bounds_violated > 0) {
    warning(sprintf(
      "cell-mode bounds: %d weights fell outside [%.3f, %.3f] due to skewed base weights within cells. Consider bounds_mode='unit' for strict per-observation bounds.",
      calib_result$n_bounds_violated, min_weight, max_weight))
  }
  if (!is.null(calib_result$n_bounds_clamped) && calib_result$n_bounds_clamped > 0) {
    warning(sprintf(
      "unit-mode bounds: %d weights clamped to [%.3f, %.3f] during per-cell water-filling.",
      calib_result$n_bounds_clamped, min_weight, max_weight))
  }

  # Enum: RK_ALG_AUTO=0, RK_ALG_IEPPA=1, RK_ALG_LBFGSB=2, RK_ALG_RAKING=3
  alg_names <- c("", "ieppa", "lbfgsb", "raking")  # index 0 (auto) removed from user API
  alg_used  <- alg_names[calib_result$algorithm_used + 1L]

  if (!attach_weights) {
    attr(weights, "result") <- calib_result
    return(weights)
  }

  col <- if (!is.null(weight_column)) weight_column else "weights"
  data[[col]] <- weights
  attr(data, "algorithm") <- alg_used
  attr(data, "result") <- calib_result
  data
}

# --- Helpers (each <= 15 lines, independently testable) ---

parse_target <- function(target, target_map) {
  if (!is.data.frame(target)) {
    if (!is.list(target))
      stop("target must be a named list of named numeric vectors or a data frame.")
    return(target)
  }
  if (!is.null(target_map)) {
    vcol <- target_map[["variable"]]; lcol <- target_map[["level"]]
    pcol <- target_map[["proportion"]]
  } else if (all(c("variable","level","proportion") %in% names(target))) {
    vcol <- "variable"; lcol <- "level"; pcol <- "proportion"
  } else if (ncol(target) == 3) {
    stop("target data frame has 3 columns but no 'variable'/'level'/'proportion' names. ",
         "Add column names or pass target_map=list(variable=..., level=..., proportion=...).")
  } else {
    stop("Cannot determine variable/level/proportion columns in target data frame.")
  }
  vars <- unique(target[[vcol]])
  lapply(stats::setNames(vars, vars), function(v) {
    sub <- target[target[[vcol]] == v, , drop = FALSE]
    stats::setNames(sub[[pcol]], sub[[lcol]])
  })
}

parse_bounds_mode <- function(x = c("cell", "unit")) {
  match.arg(x)
}

map_method <- function(method, verbose = 0) {
  method <- tolower(method)
  if (method %in% c("rake", "nrake")) {
    warning("method='", method, "' (IPF) not implemented; using L-BFGS-B")
    method <- "lbfgsb"
  } else if (method == "nr") {
    warning("method='nr' (Newton-Raphson) not implemented; using L-BFGS-B")
    method <- "lbfgsb"
  }
  match.arg(method, c("ieppa", "lbfgsb", "raking"))
}

parse_convergence <- function(convergence) {
  if (!is.null(convergence[["pct"]]))
    warning("convergence['pct'] is deprecated in leafblower; use convergence['absolute'].")
  if (!is.null(convergence[["absolute"]])) convergence[["absolute"]] else 1e-6
}

normalize_start_weights <- function(start_weights, n) {
  if (is.null(start_weights)) return(NULL)
  if (length(start_weights) == 1) {
    sw <- rep(as.double(start_weights), n)
  } else {
    if (length(start_weights) != n)
      stop("start_weights length must equal nrow(data)")
    sw <- as.double(start_weights)
  }
  if (sum(sw) < 1e-15)
    stop("start_weights must sum to a positive value")
  sw * length(sw) / sum(sw)
}
