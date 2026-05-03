#' Generate calibrated weights (drop-in for autumn::harvest)
#'
#' @param data A data frame containing columns to calibrate.
#' @param target A named list of named numeric vectors (variable -> proportions).
#' @param min_weight Lower bound on weights. Default 0 (no lower bound).
#' @param max_weight Upper bound on weights. Default 5.
#' @param capacity_penalty Numeric, controls how strongly capacity bounds are
#'   enforced during ALM optimization in \code{method="ieppa_soft"}. Use
#'   \code{NULL} (default) for auto-computed value (\code{M_cell / n}) which gives a
#'   balanced blend between unconstrained KL minimization and hard-clamp projection.
#'   Larger values force tighter constraint adherence; smaller values allow more
#'   temporary bound violation before the final projection enforces exact bounds.
#'   Tuning: if \code{attr(result, "result")$alm_capacity_mu_final / capacity_penalty >= 1000},
#'   the adaptive schedule hit its ceiling — increase \code{capacity_penalty} by 10x.
#'   Ignored for methods other than \code{"ieppa_soft"}.
#' @param method Calibration method. One of \code{"auto"} (default: iEPPA or
#'   raking based on M_cell/n ratio), \code{"ieppa"} (paper-faithful iEPPA),
#'   \code{"ieppa_soft"} (iEPPA with augmented Lagrangian soft capacity enforcement;
#'   better than \code{"ieppa"} on tight-bounds problems where cells hit
#'   \code{max_weight}),
#'   \code{"raking"} (IPF + water-filling box projection (KL projection, Csiszar-Tusnady 1984)), \code{"lbfgsb"}
#'   (L-BFGS-B on concave dual), \code{"sinkhorn"} (KL Bregman Dykstra),
#'   \code{"greg"} (Newton QP, Deville-Sarndal 1992), \code{"chebyshev"}
#'   (L-infinity LP via IPM), \code{"greenkhorn"} (greedy coordinate-descent IPF — Altschuler-Weed-Rigollet 2017;
#'   picks the single hardest margin per step; supports \code{accelerate=TRUE} for
#'   SQUAREM round-level acceleration. Distinct from \code{scheduler="greedy"} on raking,
#'   which sorts all K margins but still sweeps all per round),
#'   \code{"logit"} (Deville-Sarndal 1992 logit-distance Newton calibration; bounds
#'   enforced analytically — no clamping; typically converges in 10-20 Newton steps;
#'   equivalent to \code{autumn::calibrate()}).
#' @param verbose Integer verbosity: 0=silent, 1=progress, 2=debug.
#' @param max_iterations Maximum inner BCD iterations per outer step. Default 500.
#' @param start_weights Starting weights vector or NULL (uniform).
#' @param attach_weights If TRUE, return data frame with weights column. Default TRUE.
#' @param weight_column Name of weight column. Default "weights".
#' @param convergence Named list controlling the stopping criterion. Accepted keys:
#'   \itemize{
#'     \item \code{metric}: quality metric to monitor. One of \code{"max_err"}
#'       (default), \code{"mean_err"}, \code{"kl"}, \code{"chi2"},
#'       \code{"grake_norm"}, \code{"l1_weight"}, \code{"marginal_kl"}.
#'     \item \code{rule}: stopping rule. One of \code{"improvement"} (default) —
#'       stop when metric improves by less than \code{tol} relative; \code{"threshold"} —
#'       stop when metric falls below \code{tol} absolutely; \code{"plateau"} —
#'       stop when metric changes less than \code{tol} over a window.
#'     \item \code{tol}: tolerance value (default \code{0.001}).
#'     \item \code{stop_when}: \code{"any"} (default) or \code{"all"}.
#'   }
#'   Shorthand keys (each sets metric+rule+tol simultaneously):
#'   \itemize{
#'     \item \code{pct}: sets \code{metric="l1_weight"}, \code{rule="plateau"},
#'       \code{tol=pct}. Example: \code{list(pct = 0.001)}.
#'     \item \code{absolute}: sets \code{metric="max_err"}, \code{rule="threshold"},
#'       \code{tol=absolute}. Example: \code{list(absolute = 1e-4)}.
#'     \item \code{improvement}: sets \code{metric="max_err"}, \code{rule="improvement"},
#'       \code{tol=improvement}. Example: \code{list(improvement = 0.001)}.
#'   }
#'   Shorthands and explicit keys may be combined; explicit keys override shorthand
#'   defaults. \code{convergence = list()} uses the default (max_err + improvement +
#'   tol = 0.001). For \code{method="ieppa"}, the default metric is \code{marginal_kl}
#'   (calibration quality: Σ_k Σ_j t_kj log(t_kj/achieved_kj)). For other methods: \code{max_err}.
#'
#'   \strong{chi2 cross-solver note:} chi2 is not directly comparable
#'   across methods. iEPPA uses unnormalized cell mass as \code{W_total};
#'   raking and lbfgsb use \code{n}. Use chi2 as a convergence criterion
#'   within one method; do not compare values across methods.
#' @param sor Named list for SOR adaptive under-relaxation (iEPPA and raking).
#'   \code{NULL} disables SOR. Keys:
#'   \itemize{
#'     \item \code{auto}: logical, default \code{TRUE}.
#'     \item \code{omega_min}: lower bound on relaxation factor, default \code{0.3}.
#'     \item \code{omega}: fixed relaxation factor (disables auto-adapt).
#'     \item \code{burnin}: iterations before adaptation starts, default \code{20}.
#'   }
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
#' @param jacobi_sweep Logical. Use Jacobi (frozen cell_lf snapshot) instead of
#'   Gauss-Seidel (incremental cell_lf updates) for the iEPPA log-path margin sweep.
#'   Default \code{FALSE} (Gauss-Seidel). May improve wall-time at very large
#'   \code{M_cell} (>1M unique cells) where scattered GS writes exceed last-level
#'   cache capacity. No benefit observed at M_cell < 100K. Ignored for linear-path
#'   problems (\code{n/M_cell < 2}).
#' @param accelerate Logical. Enable Safeguarded Regularized Anderson Acceleration
#'   (SRAA-m, window m=5) for \code{method="raking"}, \code{"greenkhorn"},
#'   \code{"ieppa"}, and \code{"ieppa_soft"}. Default \code{FALSE}.
#'
#'   \strong{For ieppa and ieppa_soft}, SRAA-m operates in log-factor space
#'   (dimension n_cats_total ≈ 50–500, not M_cell), so history matrices are
#'   small. Expected benefit: ≥30\% fewer outer iterations on tight-bounds
#'   problems (K≥6, max_weight<3, or skewed margins).
#'
#'   \strong{Behavioral changes when accelerate=TRUE for ieppa/ieppa_soft:}
#'   \itemize{
#'     \item \strong{SOR disabled}: adaptive under-relaxation (SOR) is turned
#'       off; omega stays fixed at omega_init (default 1.0). Use
#'       \code{accelerate=FALSE} if adaptive SOR is required.
#'     \item \strong{Greedy scheduler downgraded}: if \code{scheduler="greedy"},
#'       it is silently changed to \code{"round_robin"} (logged at
#'       \code{verbose >= 1}). AA requires a deterministic sweep map.
#'     \item \strong{res$iterations}: reports F-evals consumed, not BCD outer
#'       iterations (AA-accepted step = 2 F-evals; plain step = 1).
#'     \item \strong{res$aa_accepted_count}: number of AA steps accepted this
#'       solve (0 when accelerate=FALSE). Access via
#'       \code{attr(result, "result")$aa_accepted_count}.
#'     \item SRAA history is reset at each homotopy level boundary and on
#'       linear-to-log-path fallback.
#'   }
#' @param add_na_proportion Not supported in v1; raises error if TRUE.
#' @param auto_collapse Not supported in v1; raises error if TRUE.
#' @param collapse_vars Not supported in v1; raises error if TRUE.
#' @param target_map Passed through for data-frame target format handling.
#' @param design_weights Optional design weights vector. When non-NULL and
#'   \code{start_weights} is NULL, used as starting weights (normalized to
#'   mean=1). Length must equal \code{nrow(data)}.
#' @param ... Additional arguments ignored.
#' @return data frame with weights column if \code{attach_weights=TRUE}, else a numeric
#'   vector of length \code{n}. In both cases the object carries attributes:
#'   \describe{
#'     \item{\code{result}}{Named list of solver diagnostics. Key fields:
#'       \itemize{
#'         \item \code{status}: integer status code. 0=converged (RK_OK);
#'           1=did not converge (RK_ERR_NOCONV); 2=infeasible (RK_ERR_INFEAS);
#'           3=bad argument (RK_ERR_BADARG); 4=budget exhausted — loss still
#'           decreasing, increase \code{max_iterations} (RK_ERR_BUDGET);
#'           5=loss plateau at constrained optimum, weights are valid (RK_ERR_STALL).
#'         \item \code{iterations}: number of outer iterations completed.
#'         \item \code{max_error}: maximum marginal error at the returned weights.
#'         \item \code{l1_weight_change}: L1-normalized weight change Σ|Δw|/n.
#'           For \code{method="lbfgsb"} this is start-to-final (batch solver);
#'           for iEPPA and raking it is iteration-to-iteration.
#'         \item \code{grake_norm}: survey-grake normalized residual
#'           max_k|misfit|/(1+|pop|).
#'         \item \code{convergence_used}: Named list with \code{metric}, \code{rule},
#'           \code{tol}, and \code{fired_at_iter} documenting which criterion fired.
#'         \item \code{best_error}: minimum marginal error seen across all iterates.
#'         \item \code{best_weights}: numeric vector (length \code{n}, sum normalized to
#'           \code{n}) at the iterate with minimum observed marginal error. When
#'           \code{best_error} is \code{Inf} (solver exited before the first convergence
#'           check), \code{best_weights} is all-zero. Guard:
#'           \code{if (is.finite(attr(r, "result")$best_error))} before use.
#'         \item \code{convergence_used$convergence_reason}: Character.
#'           Why the solver exited: \code{"criterion"} (improvement criterion satisfied),
#'           \code{"budget"} (budget exhausted — increase max_iterations),
#'           \code{"stall_kl"} (weight KL plateau — at constrained KL minimum),
#'           \code{"stall_wchange"} (SRAA-m weight-change plateau — at constrained optimum),
#'           \code{"infeasible"}, \code{"error"}, or \code{"legacy"}.
#'         \item \code{alm_capacity_mu_final}: final ALM penalty after adaptive scaling (\code{0} if not \code{ieppa_soft}).
#'         \item \code{alm_n_growth_events}: adaptive penalty growth fire count.
#'         \item \code{alm_max_dual_norm}: max absolute Lagrange dual at solver exit.
#'         \item \code{alm_sum_drift}: \code{|sum(weights) - n|} after final projection (bounded by \code{1e-6 * n}).
#'       }
#'     }
#'     \item{\code{algorithm}}{Character name of the solver used.}
#'     \item{\code{iterations}}{Convenience alias for \code{result$iterations}.}
#'   }
#' @details
#' \strong{When to use \code{ieppa_soft} vs \code{ieppa}}: Use
#' \code{method="ieppa_soft"} when \code{ieppa} gives \code{max_error > 1e-3}
#' and many observations are near \code{max_weight}. \code{ieppa_soft} is
#' roughly 10-30\% slower than \code{ieppa} but finds a better constrained
#' optimum by temporarily relaxing bounds during optimization, then projecting
#' to exact feasibility at exit.
#'
#' \strong{Greenkhorn}: \code{method="greenkhorn"} implements the
#' Altschuler-Weed-Rigollet (2017) greedy coordinate-descent Sinkhorn variant.
#' At each step \code{greenkhorn} selects the single margin with the largest
#' violation and updates only that row/column — unlike standard raking which
#' sweeps all K margins every round. Pass \code{accelerate=TRUE} to enable
#' SRAA-m outer-loop acceleration on top of greenkhorn.
#'
#' \strong{Logit calibration}: \code{method="logit"} implements Deville-Sarndal
#' (1992) logit-distance Newton calibration. Bounds are enforced analytically
#' through the logit link — no post-hoc weight clamping is required. The logit
#' method typically converges in 10–20 Newton steps and is equivalent to
#' \code{autumn::calibrate()} with a logit distance function.
#' @examples
#' \dontrun{
#' df  <- data.frame(sex = factor(sample(c("M","F"), 500, TRUE)))
#' tgt <- list(sex = c(M = 0.5, F = 0.5))
#'
#' # Greenkhorn (greedy coordinate-descent IPF)
#' r_grk <- harvest(df, tgt, method = "greenkhorn")
#'
#' # Logit-distance Newton calibration (Deville-Sarndal 1992)
#' r_logit <- harvest(df, tgt, method = "logit")
#' }
#' @export
harvest <- function(
  data,
  target,
  min_weight       = 0,
  max_weight       = 5,
  capacity_penalty = NULL,
  alm_penalty      = NULL,
  method           = "ieppa",
  verbose          = 0,
  max_iterations   = 500,
  start_weights    = NULL,
  attach_weights   = TRUE,
  weight_column    = "weights",
  convergence      = list(),
  sor              = list(auto = TRUE, omega_min = 0.3),
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
  jacobi_sweep     = FALSE,
  newton_tsvd_ratio = 1e-8,
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
  conv    <- parse_convergence(convergence)
  # iEPPA only: marginal_kl is the calibration-quality loss (best_iter criterion).
  # Raking and other solvers keep their own defaults — raking minimizes weight KL
  # via IPF, not marginal KL.
  if (method %in% c("ieppa", "ieppa_soft") &&
      is.null(convergence[["metric"]]) &&
      is.null(convergence[["criterion"]]) &&
      is.null(convergence[["improvement"]]) &&
      is.null(convergence[["pct"]]) &&
      is.null(convergence[["absolute"]])) {
    conv$metric <- "marginal_kl"
  }
  # Sinkhorn minimizes weight KL — override default metric when user hasn't specified one.
  if (method == "sinkhorn" &&
      is.null(convergence[["metric"]]) &&
      is.null(convergence[["criterion"]]) &&
      is.null(convergence[["improvement"]]) &&
      is.null(convergence[["pct"]]) &&
      is.null(convergence[["absolute"]])) {
    conv$metric <- "kl"
  }
  if (!is.null(capacity_penalty)) {
    if (!is.numeric(capacity_penalty) || length(capacity_penalty) != 1L ||
        !is.finite(capacity_penalty) || capacity_penalty <= 0) {
      stop("capacity_penalty must be NULL (auto) or a positive finite scalar; got: ",
           deparse(capacity_penalty), call. = FALSE)
    }
    if (capacity_penalty > 1e15) {
      stop("capacity_penalty must be NULL (auto) or a positive finite scalar; got: ",
           deparse(capacity_penalty), call. = FALSE)
    }
    if (capacity_penalty < 1e-15) {
      warning("capacity_penalty=", capacity_penalty,
              " is below recommended range; constraint enforcement may be ineffective",
              call. = FALSE)
    }
  }

  if (!is.null(alm_penalty)) {
    if (!is.numeric(alm_penalty) || length(alm_penalty) != 1L ||
        !is.finite(alm_penalty) || alm_penalty <= 0) {
      stop("alm_penalty must be NULL (disabled) or a positive finite scalar; got: ",
           deparse(alm_penalty), call. = FALSE)
    }
    if (alm_penalty > 1e15) {
      stop("alm_penalty must be NULL (disabled) or a positive finite scalar; got: ",
           deparse(alm_penalty), call. = FALSE)
    }
    if (alm_penalty < 1e-15) {
      warning("alm_penalty=", alm_penalty,
              " is below recommended range; objective penalty may be ineffective",
              call. = FALSE)
    }
  }

  # Epic-H WH-e: newton_tsvd_ratio validation (newton_kl only; ignored elsewhere by C side).
  if (!is.numeric(newton_tsvd_ratio) || length(newton_tsvd_ratio) != 1L ||
      !is.finite(newton_tsvd_ratio) || newton_tsvd_ratio <= 0) {
    stop("newton_tsvd_ratio must be a positive finite scalar; got: ",
         deparse(newton_tsvd_ratio), call. = FALSE)
  }

  sor_cfg <- parse_sor(sor)
  if (isTRUE(accelerate) && !method %in% c("raking", "greenkhorn", "ieppa", "ieppa_soft"))
    warning("accelerate=TRUE is only supported for method='raking', 'greenkhorn', 'ieppa', or 'ieppa_soft'; ignoring for method='",
            method, "'", call. = FALSE)
  accelerate_bool <- isTRUE(accelerate) && method %in% c("raking", "greenkhorn", "ieppa", "ieppa_soft")

  # design_weights: used as start_weights when supplied (normalized to mean=1 by normalize_start_weights)
  if (!is.null(design_weights)) {
    if (!is.null(start_weights))
      warning("leafblower: both design_weights and start_weights supplied; design_weights ignored")
    else
      start_weights <- design_weights
  }
  sw_vec  <- normalize_start_weights(start_weights, nrow(data))

  bounds_mode_char <- parse_bounds_mode(bounds_mode)
  bounds_mode_int  <- match(bounds_mode_char, c("cell", "unit")) - 1L

  # Overlay arg resolution
  scheduler    <- match.arg(scheduler)
  eta_schedule <- match.arg(eta_schedule)
  # R8: raking-only guard. Greedy re-sorts within each F_eval for raking (non-stationary).
  # greenkhorn sorts ONCE at F_eval entry (stationary); scheduler param irrelevant to it.
  if (accelerate_bool && method == "raking" && scheduler == "greedy")
    scheduler <- "round_robin"

  if (!is.null(capacity_penalty) && !grepl("ieppa_soft", method, fixed = TRUE)) {
    warning("capacity_penalty is only used by method='ieppa_soft'; ignored for method='",
            method, "'", call. = FALSE)
  }

  # Ignored-param verbose notes
  # enforce_mean is always TRUE: normalization is unconditional (line ~86).
  ignored <- c("select_params", "select_function", "error_function",
                "adaptive_order", "enforce_mean")
  supplied_ignored <- intersect(ignored, names(match.call(expand.dots = FALSE)))
  if (verbose >= 2 && length(supplied_ignored) > 0)
    message("[leafblower] Ignoring autumn params: ", paste(supplied_ignored, collapse = ", "))

  # CalibMetric: 0=MAX_ERR 1=MEAN_ERR 2=KL 3=CHI2 4=GRAKE_NORM 5=L1_WEIGHT 6=MARGINAL_KL
  metric_int    <- c(max_err = 0L, mean_err = 1L, kl = 2L, chi2 = 3L,
                     grake_norm = 4L, l1_weight = 5L,
                     marginal_kl = 6L,
                     pct = 5L)  # pct is legacy alias for l1_weight
  # CalibRule: 0=THRESHOLD 1=IMPROVEMENT 2=PLATEAU
  rule_int      <- c(threshold = 0L, improvement = 1L, plateau = 2L)
  stop_when_int <- c(any = 0L, all = 1L)

  raw <- .Call("C_rk_calibrate",
               data,
               target,
               as.double(min_weight),
               as.double(max_weight),
               as.character(method),
               as.integer(verbose),
               as.integer(max_iterations),
               sw_vec,
               if (is.null(capacity_penalty)) -1.0 else as.double(capacity_penalty),  # 9: capacity_penalty
               if (is.null(alm_penalty)) -1.0 else as.double(alm_penalty),  # 10: alm_penalty
               as.double(if (conv$absolute_tol > 0) conv$absolute_tol else 1e-6),  # slot 11: legacy tol_abs
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
               ## Convergence config (WU-A)
               as.double(conv$pct_tol),
               as.double(conv$absolute_tol),
               as.integer(metric_int[[conv$metric]]),
               as.integer(rule_int[[conv$rule]]),
               as.integer(stop_when_int[[conv$stop_when]]),
               ## SOR config (WU-A)
               as.integer(sor_cfg$enabled),
               as.integer(sor_cfg$auto),
               as.double(sor_cfg$omega_init),
               as.double(sor_cfg$omega_min),
               as.double(sor_cfg$omega_fixed),
               as.integer(sor_cfg$burnin),
               ## SRAA-m accelerate flag
               as.integer(accelerate_bool),
               ## Jacobi log-path sweep flag
               as.integer(isTRUE(jacobi_sweep)),
               ## Epic-H WH-e: newton_kl TSVD truncation ratio (default 1e-8)
               as.double(newton_tsvd_ratio),
               PACKAGE = "leafblower")

  weights <- raw$weights
  calib_result    <- raw$result

  # WU-D: nest SOR diagnostics under $sor for clean namespace.
  # The C bridge always returns sor_min_omega and sor_n_damped as flat fields;
  # wrap them here so callers use result$sor$min_omega and result$sor$n_damped.
  calib_result$sor <- list(
    min_omega = calib_result$sor_min_omega,
    n_damped  = calib_result$sor_n_damped
  )
  calib_result$sor_min_omega <- NULL
  calib_result$sor_n_damped  <- NULL

  # WU-E2: nest convergence diagnostics under $convergence_used for clean namespace.
  # metric_names and rule_names mirror CalibMetric/CalibRule enum order in leafblower.h.
  .metric_names <- c("max_err", "mean_err", "kl", "chi2", "grake_norm", "l1_weight", "marginal_kl")
  .rule_names   <- c("threshold", "improvement", "plateau")
  # C1: guard +1L indexing — integer index from C may be NA or out of range.
  .safe_lookup <- function(v, idx) {
    if (is.integer(idx) && !is.na(idx) && idx >= 0L && idx < length(v)) v[idx + 1L]
    else NA_character_
  }
  calib_result$convergence_used <- list(
    metric           = .safe_lookup(.metric_names, calib_result$convergence_metric),
    rule             = .safe_lookup(.rule_names,   calib_result$convergence_rule),
    tol              = calib_result$convergence_tol,
    fired_at_iter    = calib_result$convergence_iter,
    solver_objective = calib_result$solver_objective,   # Task 2: solver mathematical objective
    minimized_metric = .safe_lookup(.metric_names, calib_result$convergence_minimized_metric),
    convergence_reason = {
      s <- calib_result$status
      if      (is.null(s) || is.na(s))              NA_character_
      else if (s == 0L)                              "criterion"
      else if (s == 4L)                              "budget"
      else if (s == 5L && isTRUE(accelerate_bool))  "stall_wchange"
      else if (s == 5L)                              "stall_kl"
      else if (s == 2L)                              "infeasible"
      else if (s == 3L)                              "error"
      else                                           "legacy"
    }
  )
  calib_result$convergence_metric           <- NULL
  calib_result$convergence_rule             <- NULL
  calib_result$convergence_tol              <- NULL
  calib_result$convergence_iter             <- NULL
  calib_result$solver_objective             <- NULL
  calib_result$convergence_minimized_metric <- NULL

  # Check hard-stop statuses before normalization: status 2/3 mean weights are
  # meaningless; normalizing near-zero weights before stopping produces NaN.
  if (calib_result$status == 2L)
    stop("leafblower: ", if (nchar(calib_result$message) > 0) calib_result$message
         else "infeasible problem")
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

  if (calib_result$status == 4L)
    warning("leafblower: budget exhausted — weights reflect best iterate; ",
            "increase max_iterations if further improvement is needed")
  if (calib_result$status == 5L && isTRUE(accelerate_bool))
    warning("leafblower: SRAA-m weight-change plateau — at constrained optimum; ",
            "weights are valid; no further improvement is achievable")
  if (calib_result$status == 5L && !isTRUE(accelerate_bool))
    warning("leafblower: loss function plateau — at constrained optimum given bounds; ",
            "weights are valid; no further improvement is achievable")
  if (calib_result$status == 1L)
    warning("leafblower: did not converge (legacy status code from solver not yet updated)")

  # Stall detection: PCT converged (status=0) but max_error >> pct_tol
  # signals infeasible problem. Threshold: 10x pct_tol.
  # Threshold derivation: well-posed problems have errRp/pct_change ratio 1-5x;
  # infeasible stalls show 100x+; 10x cleanly separates the two regimes.
  if (calib_result$status == 0L &&
      conv$metric %in% c("max_err", "mean_err") &&
      !is.null(conv$pct_tol) && conv$pct_tol > 0 &&
      !is.null(calib_result$max_error) &&
      calib_result$max_error > 10 * conv$pct_tol) {
    warning(sprintf(
      "leafblower: PCT convergence stall: max_error=%.3g >> pct_tol=%.3g (%.0fx). Possible contradictory or infeasible targets.",
      calib_result$max_error, conv$pct_tol,
      calib_result$max_error / conv$pct_tol),
      call. = FALSE)
  }

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

  # algorithm_used is emitted as a character string by C++ (single source of truth).
  alg_used <- calib_result$algorithm_used

  # Quality-check warning: greg may be unreliable when max_err > 5%
  if (alg_used == "greg" &&
      !is.null(calib_result$max_error) &&
      is.finite(calib_result$max_error) &&
      calib_result$max_error > 0.05) {
    warning(sprintf(
      paste0("greg converged but max_err=%.4g (>5%%). ",
             "greg may be unreliable for K=%d margins or tight bounds ",
             "(max_weight=%.4g). Consider method='raking' or 'ieppa'."),
      calib_result$max_error, length(target), max_weight),
      call. = FALSE)
  }

  # Tier-1/2 calibration quality metrics.
  # Computed post-solver at O(nK) cost.
  # Tier 1: margin_kl = KL(T_kj || W_kj) summed over all margins k and levels j,
  #         where W_kj = sum of calibrated weights on category j / total weight.
  # Tier 2: weight_kl  = KL(w/Z || 1/n) = mean(wn * log(wn)) where wn = n*w/Z,
  #         design_effect   = n * sum(w^2) / sum(w)^2 (DEFF),
  #         effective_observations = sum(w)^2 / sum(w^2) (Kish n_eff).
  if (length(weights) > 0L && length(target) > 0L &&
      is.finite(sum(weights)) && sum(weights) > 0) {
    Z  <- sum(weights)
    n  <- length(weights)
    wn <- weights * n / Z                         # normalize to mean=1
    deff <- n * sum(weights^2) / Z^2
    calib_result$design_effect          <- deff
    calib_result$effective_observations <- if (is.finite(deff) && deff > 0) n / deff else NA_real_
    calib_result$weight_kl              <- sum(wn * log(pmax(wn, 1e-15))) / n
    calib_result$margin_kl              <- tryCatch(
      sum(sapply(names(target), function(k) {
        T_k   <- target[[k]]
        obs_k <- data[[k]]
        valid <- !is.na(obs_k)                    # exclude NA-coded observations
        w_v   <- weights[valid]
        Z_k   <- sum(w_v)
        if (Z_k == 0) return(NA_real_)
        W_k <- tapply(w_v, droplevels(obs_k[valid]), sum) / Z_k
        # categories in T_k with T_k > 0 but absent from W_k → infeasible margin
        if (any(T_k[setdiff(names(T_k), names(W_k))] > 0)) return(Inf)
        common <- intersect(names(T_k), names(W_k))
        T_sub  <- T_k[common]; W_sub <- W_k[common]
        pos    <- T_sub > 0                       # limit: 0 * log(0/W) = 0
        if (!any(pos)) return(0)
        sum(T_sub[pos] * log(T_sub[pos] / pmax(W_sub[pos], 1e-15)))
      }), na.rm = FALSE),
      error = function(e) NA_real_
    )
  } else {
    calib_result$design_effect          <- NA_real_
    calib_result$effective_observations <- NA_real_
    calib_result$weight_kl              <- NA_real_
    calib_result$margin_kl              <- NA_real_
  }

  if (!attach_weights) {
    attr(weights, "result")     <- calib_result
    attr(weights, "algorithm")  <- alg_used
    attr(weights, "iterations") <- calib_result$iterations
    return(weights)
  }

  col <- if (!is.null(weight_column)) weight_column else "weights"
  data[[col]] <- weights
  attr(data, "algorithm") <- alg_used
  attr(data, "result") <- calib_result
  attr(data, "iterations") <- calib_result$iterations
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
  match.arg(method, c("auto", "ieppa", "ieppa_soft", "lbfgsb", "raking", "sinkhorn", "chebyshev", "greg", "greenkhorn", "logit", "newton_kl"))
}

parse_convergence <- function(convergence) {
  if (!is.null(convergence) && !is.list(convergence))
    stop("convergence must be a named list or NULL (e.g. list(pct = 0.001))")
  valid_keys <- c("pct", "absolute", "metric", "criterion", "rule", "stop_when",
                  "tol", "improvement")
  bad <- setdiff(names(convergence), valid_keys)
  if (length(bad))
    stop(sprintf("Unknown convergence key(s): %s. Valid keys: %s",
                 paste(bad, collapse = ", "),
                 paste(valid_keys, collapse = ", ")))
  `%||%` <- function(a, b) if (is.null(a)) b else a
  explicit_pct         <- !is.null(convergence[["pct"]])
  explicit_abs         <- !is.null(convergence[["absolute"]])
  explicit_improvement <- !is.null(convergence[["improvement"]])
  explicit_tol         <- !is.null(convergence[["tol"]])

  # C2: improvement= and absolute= together without stop_when is ambiguous — error.
  if (explicit_improvement && explicit_abs && is.null(convergence[["stop_when"]])) {
    stop("convergence: 'improvement' and 'absolute' cannot be combined without ",
         "'stop_when'. Use stop_when = 'any' or 'all' to fire on either or both.")
  }

  # Shorthand: improvement=X -> max_err + improvement rule + X as pct_tol
  if (explicit_improvement) {
    tol_val <- convergence[["improvement"]]
    return(list(pct_tol = as.double(tol_val), absolute_tol = 0.0,
                metric = "max_err", rule = "improvement", stop_when = "any"))
  }

  pct_tol <- if (explicit_pct) convergence[["pct"]]
             else if (!explicit_abs) 0.001
             else 0.0
  absolute_tol <- convergence[["absolute"]] %||% 0.0

  # "criterion" is a legacy alias for "metric" (backward compat)
  metric_raw <- convergence[["metric"]] %||% convergence[["criterion"]] %||%
                (if (explicit_pct) "pct" else "max_err")
  metric <- match.arg(metric_raw,
    c("max_err", "mean_err", "kl", "chi2", "grake_norm", "l1_weight", "marginal_kl", "pct"))

  # pct is autumn/anesrake compatible: stops when Σ|Δw| STOPS IMPROVING (plateau).
  # Default its rule to "plateau" when pct is specified without an explicit rule.
  # absolute= with no explicit rule maps to "threshold" (hard stopping criterion).
  rule_explicit <- !is.null(convergence[["rule"]])
  rule_default  <- if (!rule_explicit && explicit_pct) "plateau"
                   else if (!rule_explicit && explicit_abs && !explicit_pct) "threshold"
                   else "improvement"
  rule_raw      <- convergence[["rule"]] %||% rule_default
  rule         <- match.arg(rule_raw, c("threshold", "improvement", "plateau"))

  # "tol" shorthand: overrides pct_tol for threshold rule, pct_tol otherwise.
  if (explicit_tol) {
    tol_val <- as.double(convergence[["tol"]])
    if (rule == "threshold") {
      absolute_tol <- tol_val
      pct_tol      <- 0.0
    } else {
      pct_tol <- tol_val
    }
    if (rule == "plateau" && (tol_val <= 0 || tol_val >= 1))
      stop("convergence$tol must be in (0,1) for rule='plateau'")
  }

  stop_when <- match.arg(convergence[["stop_when"]] %||% "any", c("any", "all"))
  list(pct_tol = pct_tol, absolute_tol = absolute_tol,
       metric = metric, rule = rule, stop_when = stop_when)
}

parse_sor <- function(sor) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  if (is.null(sor)) {
    return(list(enabled = 0L, auto = 0L, omega_init = 1.0,
                omega_min = 0.3, omega_fixed = -1.0, burnin = 20L))
  }
  valid_keys <- c("auto", "omega_min", "omega", "omega_init", "burnin")
  bad <- setdiff(names(sor), valid_keys)
  if (length(bad))
    stop(sprintf("Unknown sor key(s): %s. Valid keys: %s",
                 paste(bad, collapse = ", "),
                 paste(valid_keys, collapse = ", ")))
  list(
    enabled     = 1L,
    auto        = if (isTRUE(sor[["auto"]])) 1L else 0L,
    omega_init  = as.double(sor[["omega_init"]] %||% 1.0),
    omega_min   = as.double(sor[["omega_min"]] %||% 0.3),
    omega_fixed = as.double(sor[["omega"]] %||% -1.0),
    burnin      = as.integer(sor[["burnin"]] %||% 20L)
  )
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
