#' Generate calibrated weights (drop-in for autumn::harvest)
#'
#' @param data A data frame containing columns to calibrate.
#' @param target A named list of named numeric vectors (variable -> proportions).
#' @param min_weight Lower bound on weights. Default 0 (no lower bound).
#' @param max_weight Upper bound on weights. Default 5.
#' @param capacity_penalty Numeric, controls how strongly capacity bounds are
#'   enforced during ALM optimization in \code{method="oris_soft"}. Use
#'   \code{NULL} (default) for auto-computed value (\code{estimate\_M\_cell() / n},
#'   capped at 1.0 for K>8) which gives a balanced blend between unconstrained KL
#'   minimization and hard-clamp projection.
#'   Larger values force tighter constraint adherence; smaller values allow more
#'   temporary bound violation before the final projection enforces exact bounds.
#'   Tuning: if \code{attr(result, "result")$alm_capacity_mu_final / capacity_penalty >= 1000},
#'   the adaptive schedule hit its ceiling — increase \code{capacity_penalty} by 10x.
#'   Ignored for methods other than \code{"oris_soft"}.
#' @param method Calibration method (ORIS: Over-Relaxed Iterative Scaling). One of
#'   \code{"auto"} (default: ORIS or raking based on M_cell/n ratio),
#'   \code{"oris"} (paper-faithful ORIS), \code{"oris_soft"} (ORIS with augmented
#'   Lagrangian soft capacity enforcement; better than \code{"oris"} on tight-bounds
#'   problems where cells hit \code{max_weight}),
#'   \code{"raking"} (IPF + water-filling box projection (KL projection, Csiszar-Tusnady 1984)), \code{"sinkhorn"} (KL Bregman Dykstra),
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
#'   tol = 0.001). For \code{method="oris"}, the default metric is \code{marginal_kl}
#'   (calibration quality: Σ_k Σ_j t_kj log(t_kj/achieved_kj)). For other methods: \code{max_err}.
#'
#'   \strong{chi2 cross-solver note:} chi2 is not directly comparable
#'   across methods. ORIS uses unnormalized cell mass as \code{W_total};
#'   raking uses \code{n}. Use chi2 as a convergence criterion
#'   within one method; do not compare values across methods.
#'
#'   \strong{chebyshev note:} \code{method="chebyshev"} ignores \code{rule}
#'   (and \code{metric}/\code{stop_when}). The interior-point solver stops on its
#'   own complementarity-gap criterion, falling back to the \code{pct}/\code{absolute}
#'   tolerance on the max marginal error. Other \code{convergence} keys are no-ops
#'   for this method.
#' @param sor Named list for SOR adaptive under-relaxation (ORIS and raking).
#'   Default \code{NULL} disables SOR. Keys:
#'   \itemize{
#'     \item \code{auto}: logical, default \code{TRUE}.
#'     \item \code{omega_min}: lower bound on relaxation factor, default \code{0.3}.
#'     \item \code{omega_max}: upper bound (recovery ceiling) on relaxation factor,
#'       default \code{1.5}. Values in (1, 2) enable over-relaxation, which reduces
#'       iteration count; global convergence is guaranteed for all omega in (0, 2)
#'       (Thibault 2021). ORIS only; ignored for greedy scheduler.
#'     \item \code{omega}: fixed relaxation factor (disables auto-adapt).
#'     \item \code{burnin}: iterations before adaptation starts, default \code{20}.
#'     \item \code{omega_mode_id}: omega adaptation strategy. \code{0} = heuristic
#'       (0.7 damp / 1.05 grow), \code{1} = fixed (jump to \code{omega_max}),
#'       \code{2} = iterate-change (free-coordinate \eqn{\|\Delta X_{\rm free}\|^2}
#'       estimator, \strong{default}; feasibility-agnostic; e18t.9 SHIP — 240 vs 350
#'       iters on T2 unconstrained, 50 vs 140 on bounded stepstone).
#'       String aliases \code{"heuristic"}, \code{"fixed"}, \code{"spectral"} are accepted.
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
#' @param accelerate Logical. Enable Safeguarded Regularized Anderson Acceleration
#'   (SRAA-m, window m=5) for \code{method="raking"}, \code{"greenkhorn"},
#'   \code{"oris"}, and \code{"oris_soft"}. Default \code{FALSE}.
#'
#'   \strong{For oris and oris_soft}, SRAA-m operates in log-factor space
#'   (dimension n_cats_total ≈ 50–500, not M_cell), so history matrices are
#'   small. Expected benefit: ≥30\% fewer outer iterations on tight-bounds
#'   problems (K≥6, max_weight<3, or skewed margins).
#'
#'   \strong{Behavioral changes when accelerate=TRUE for oris/oris_soft:}
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
#' @param ridge_lambda Tikhonov ridge penalty on the dual lambda (default 0.0 = off).
#'   Adds \code{ridge_lambda} to the diagonal of the Hessian (newton_kl) or
#'   normal-equations matrix (greg) before factorization. Hardens near-singular
#'   systems from cell-level sparseness. Ignored by all other solvers.
#' @param add_na_proportion Logical (default FALSE). When TRUE, NA observations
#'   in each margin column are encoded as an explicit \code{"NA"} bin: existing
#'   targets are renormalized by \code{(1 - na_frac)} and an \code{"NA"} target
#'   equal to \code{na_frac} is injected, so that the NA observations receive
#'   positive calibration weights rather than being excluded. An error is raised
#'   if a margin is entirely NA (\code{na_frac == 1}).
#'   \strong{Known limitation:} a column that simultaneously contains real
#'   \code{NA} values \emph{and} a literal factor level or character value named
#'   \code{"NA"} will collide — both will be mapped to the injected NA bin.
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
#'           for ORIS and raking it is iteration-to-iteration.
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
#'         \item \code{alm_capacity_mu_final}: final ALM penalty after adaptive scaling (\code{0} if not \code{oris_soft}).
#'         \item \code{alm_n_growth_events}: adaptive penalty growth fire count.
#'         \item \code{alm_max_dual_norm}: max absolute Lagrange dual at solver exit.
#'         \item \code{alm_sum_drift}: \code{|sum(weights) - n|} after final projection (bounded by \code{1e-6 * n}).
#'         \item \code{sraa_demoted}: logical; \code{TRUE} iff SRAA-m
#'           acceleration was requested with \code{scheduler="greedy"} and the
#'           greedy scheduler was demoted to round-robin (greedy is incompatible
#'           with SRAA's fixed-point geometry). \code{FALSE} for non-oris /
#'           non-raking solvers and whenever no demotion occurred.
#'       }
#'     }
#'     \item{\code{algorithm}}{Character name of the solver used.}
#'     \item{\code{iterations}}{Convenience alias for \code{result$iterations}.}
#'   }
#' @details
#' \strong{When to use \code{oris_soft} vs \code{oris}}: Use
#' \code{method="oris_soft"} when \code{oris} gives \code{max_error > 1e-3}
#' and many observations are near \code{max_weight}. \code{oris_soft} is
#' roughly 10-30\% slower than \code{oris} but finds a better constrained
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
  method           = "oris",
  verbose          = 0,
  max_iterations   = 500,
  start_weights    = NULL,
  attach_weights   = TRUE,
  weight_column    = "weights",
  convergence      = list(),
  sor              = NULL,
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
  accelerate       = FALSE,
  add_na_proportion = FALSE,
  auto_collapse    = FALSE,
  collapse_vars    = NULL,
  target_map       = NULL,
  design_weights   = NULL,
  newton_tsvd_ratio = 1e-8,
  ridge_lambda = 0.0,
  ...
) {
  # RVAL.2: warn on unknown ... args (typos / removed params)
  dots <- list(...)
  if (length(dots) > 0L)
    warning("harvest: unknown argument(s) ignored: ",
            paste(names(dots), collapse = ", "), call. = FALSE)

  # Not-in-v1 hard stops
  target  <- parse_target(target, target_map)

  # --- 81bx: auto_collapse — merge rare categories into __other__ ---
  if (isTRUE(auto_collapse) || (!is.null(collapse_vars) && !isFALSE(auto_collapse))) {
    vars_to_collapse <- if (!is.null(collapse_vars)) collapse_vars else names(target)
    data_modified    <- FALSE
    data_local       <- data
    for (v in intersect(vars_to_collapse, names(target))) {
      rare <- names(which(
        vapply(names(target[[v]]), function(lv) {
          target[[v]][[lv]] < 0.01 ||
          sum(data_local[[v]] == lv, na.rm = TRUE) < 30L
        }, logical(1))
      ))
      rare <- setdiff(rare, "NA")   # never collapse the NA bin
      if (length(rare) == 0L) next
      other_mass <- sum(unlist(target[[v]][rare]))
      target[[v]] <- target[[v]][setdiff(names(target[[v]]), rare)]
      # target[[v]] is a named numeric vector, not a list: `[[missing]]` errors
      # (subscript out of bounds) rather than returning NULL, so guard the lookup.
      existing_other <- if ("__other__" %in% names(target[[v]]))
                          target[[v]][["__other__"]] else 0
      target[[v]][["__other__"]] <- existing_other + other_mass
      # Recode observations: rare values → "__other__"
      n0  <- sum(is.na(data_local[[v]]))
      col <- as.character(data_local[[v]])
      col[col %in% rare] <- "__other__"
      data_local[[v]] <- factor(col, levels = names(target[[v]]))
      n1  <- sum(is.na(data_local[[v]]))
      if (n1 - n0 > 0L)
        warning(sprintf(
          "auto_collapse: %d row(s) of '%s' became NA (value not in collapsed target levels)",
          n1 - n0, v
        ), call. = FALSE)
      data_modified   <- TRUE
    }
    if (data_modified) data <- data_local
  }

  # --- yaye: add_na_proportion — encode NA observations as explicit "NA" bin ---
  # Tracks which margins received a NA bin so group_ids encoding can use the
  # character path with explicit NA→"NA" fill (as.character(NA) produces
  # NA_character_, NOT the string "NA"; the fill is applied in the else-branch
  # below) rather than the factor path which maps NA codes to -1L.
  .na_margins <- character(0)
  if (isTRUE(add_na_proportion)) {
    n_local <- nrow(data)
    data_names <- names(data)
    for (v in names(target)) {
      if (!v %in% data_names) next
      na_frac <- mean(is.na(data[[v]]))
      if (na_frac == 0) next
      if (na_frac == 1)
        stop(sprintf(
          "add_na_proportion: all observations are NA for margin '%s'", v),
          call. = FALSE)
      # RVAL.1: validate target sums to 1 before rescaling
      s <- sum(unlist(target[[v]]))
      if (abs(s - 1) >= 1e-6)
        stop(sprintf(
          "add_na_proportion: target for variable '%s' must sum to 1 before rescaling (observed sum = %.8g)",
          v, s), call. = FALSE)
      # Renormalize existing targets by (1 - na_frac) then add NA bin.
      # CR-F1 (dtkn.1): must stay a named NUMERIC vector — the old
      # c(lapply(...), list("NA"=)) built a LIST, so downstream margin_kl_one
      # arithmetic hit "non-numeric argument to binary operator", was swallowed
      # by tryCatch, and every add_na_proportion=TRUE run reported margin_kl=NA.
      # unlist() preserves the level names for both list- and vector-typed targets.
      target[[v]] <- c(
        unlist(target[[v]]) * (1 - na_frac),
        "NA" = na_frac
      )
      .na_margins <- c(.na_margins, v)
    }
  }

  # --- c8w1: sparseness diagnostic (pre-solve) ---
  sparse_diag <- compute_sparseness_diag(target, data,
                                         cat_threshold = 0.01, obs_threshold = 30L,
                                         na_margins = .na_margins)
  if (length(sparse_diag) > 0) {
    n_flagged <- sum(vapply(sparse_diag, length, FUN.VALUE = integer(1)))
    warning(sprintf(
      "leafblower: %d sparse categories detected (T_kj < 1%% or n_kj < 30); see result$diagnostics$sparseness",
      n_flagged), call. = FALSE)
  }

  method  <- map_method(method, verbose)
  conv    <- parse_convergence(convergence)
  # Per-method default convergence metric: use the solver's natural objective.
  # Only applied when the user has not explicitly specified any metric or tol shorthand.
  # User-provided metric/criterion/improvement/pct/absolute always take precedence.
  .no_explicit_metric <-
    is.null(convergence[["metric"]])      &&
    is.null(convergence[["criterion"]])   &&
    is.null(convergence[["improvement"]]) &&
    is.null(convergence[["pct"]])         &&
    is.null(convergence[["absolute"]])
  if (.no_explicit_metric) {
    conv$metric <- switch(method,
      # KL minimizers — marginal_kl is monotone across full sweeps (Csiszar-Tusnady).
      "oris"        = "marginal_kl",
      "oris_soft"   = "marginal_kl",
      "auto"        = "marginal_kl",  # auto routes to oris in most cases
      # Weight-KL minimizers — kl monotone by Csiszar-Tusnady; marginal_kl not
      # computed in raking/greenkhorn need_extra gate so cannot be used.
      "raking"      = "kl",
      "greenkhorn"  = "kl",
      "sinkhorn"    = "kl",
      "newton_kl"   = "kl",
      # chi2 minimizer — use its actual objective.
      "greg"        = "chi2",
      # Remaining (chebyshev minimizes L-inf = max_err; logit has no natural KL):
      conv$metric   # keep max_err default
    )
  }
  # CR-D2 (dtkn/j7x8.2): reject NA/non-finite iteration budgets and bounds before
  # they reach the C layer as INT_MIN / NaN and return garbage weights as RK_OK.
  if (length(max_iterations) != 1L || is.na(max_iterations) ||
      !is.finite(max_iterations) || max_iterations < 1) {
    stop("max_iterations must be a single positive integer; got: ",
         deparse(max_iterations), call. = FALSE)
  }
  if (length(min_weight) != 1L || is.na(min_weight) || !is.finite(min_weight)) {
    stop("min_weight must be a single finite number; got: ",
         deparse(min_weight), call. = FALSE)
  }
  if (length(max_weight) != 1L || is.na(max_weight)) {   # +Inf allowed (unbounded)
    stop("max_weight must be a single non-NA number (Inf allowed); got: ",
         deparse(max_weight), call. = FALSE)
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
  if (isTRUE(accelerate) && !method %in% c("raking", "greenkhorn", "oris", "oris_soft"))
    warning("accelerate=TRUE is only supported for method='raking', 'greenkhorn', 'oris', or 'oris_soft'; ignoring for method='",
            method, "'", call. = FALSE)
  accelerate_bool <- isTRUE(accelerate) && method %in% c("raking", "greenkhorn", "oris", "oris_soft")

  # Validate data: must be a non-empty data.frame
  if (!is.data.frame(data) || nrow(data) == 0L)
    stop("leafblower: 'data' must be a non-empty data.frame", call. = FALSE)

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
  if (is.na(bounds_mode_int))
    stop(sprintf("bounds_mode must be one of 'cell' or 'unit'; got: %s",
                 format(bounds_mode_char)), call. = FALSE)

  # Overlay arg resolution
  scheduler    <- match.arg(scheduler)
  eta_schedule <- match.arg(eta_schedule)
  # R8: raking-only guard. Greedy re-sorts within each F_eval for raking (non-stationary).
  # greenkhorn sorts ONCE at F_eval entry (stationary); scheduler param irrelevant to it.
  if (accelerate_bool && method == "raking" && scheduler == "greedy")
    scheduler <- "round_robin"

  if (!is.null(capacity_penalty) && !grepl("oris_soft", method, fixed = TRUE)) {
    warning("capacity_penalty is only used by method='oris_soft'; ignored for method='",
            method, "'", call. = FALSE)
  }

  # CalibMetric: 0=MAX_ERR 1=MEAN_ERR 2=KL 3=CHI2 4=GRAKE_NORM 5=L1_WEIGHT 6=MARGINAL_KL
  # META.2: pct alias resolved to l1_weight in parse_convergence; no duplicate entry here
  metric_int    <- c(max_err = 0L, mean_err = 1L, kl = 2L, chi2 = 3L,
                     grake_norm = 4L, l1_weight = 5L,
                     marginal_kl = 6L)
  # CalibRule: 0=THRESHOLD 1=IMPROVEMENT 2=PLATEAU
  rule_int      <- c(threshold = 0L, improvement = 1L, plateau = 2L)
  stop_when_int <- c(any = 0L, all = 1L)

  # Validate all target variables exist in data
  missing_vars <- setdiff(names(target), names(data))
  if (length(missing_vars) > 0L)
    stop("leafblower: variable(s) not found in data: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)

  margins      <- names(target)
  # Encode each margin column to 0-indexed integer (NA/-1 for OOV/NA).
  # Factor path: map factor integer codes via a precomputed level→target-index
  # table (O(nlevels) string ops once, then O(n) int array indexing — no
  # as.character() allocation). ~5× faster than the character path on large n.
  # Character path: fallback for non-factor columns via match().
  group_ids_r  <- lapply(margins, function(v) {
    col <- data[[v]]
    # yaye: if this margin has a NA bin, always use character path so that
    # NA observations can be explicitly filled to "NA" (as.character(NA) yields
    # NA_character_, not the string "NA"; the fill happens in the else-branch).
    if (is.factor(col) && !v %in% .na_margins) {
      lvl_map              <- match(levels(col), names(target[[v]])) - 1L
      lvl_map[is.na(lvl_map)] <- -1L          # OOV levels -> -1
      codes                <- as.integer(col)  # 1-indexed; NA -> NA_integer_
      gids                 <- lvl_map[codes]   # NA index -> NA result
      gids[is.na(gids)]    <- -1L
      gids
    } else {
      key <- as.character(col)
      if (v %in% .na_margins) key[is.na(key)] <- "NA"
      idx              <- match(key, names(target[[v]]))
      idx[is.na(idx)]  <- 0L
      idx - 1L
    }
  })
  # RVAL.3: warn on genuine OOV observations (gid == -1, not NA-on-NA-margin)
  for (i in seq_along(margins)) {
    v    <- margins[[i]]
    gids <- group_ids_r[[i]]
    # CR-F6 (dtkn.6): plain NAs fold to gid==-1 too (as.integer(factor) -> NA -> -1L,
    # match(NA) -> 0 -> -1L). When add_na_proportion=FALSE these are ORDINARY missing
    # data, not a vocabulary problem — exclude them so the OOV warning counts only
    # GENUINE out-of-vocabulary values (present in data, absent from the target).
    n_oov <- sum(gids == -1L & !is.na(data[[v]]))
    if (n_oov > 0L)
      warning(sprintf(
        "harvest: %d out-of-vocabulary observation(s) for variable '%s' — these levels are absent from target and will not contribute to calibration",
        n_oov, v), call. = FALSE)
  }
  cat_counts_r <- vapply(target, length, integer(1L))
  targets_r    <- lapply(target, function(t) as.double(unname(t)))
  n_obs        <- nrow(data)

  raw <- .Call("C_rk_calibrate",
               group_ids_r,
               cat_counts_r,
               targets_r,
               as.integer(n_obs),
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
               as.double(sor_cfg$omega_max),
               as.double(sor_cfg$omega_fixed),
               as.integer(sor_cfg$burnin),
               as.integer(sor_cfg$omega_mode_id),
               ## SRAA-m accelerate flag
               as.integer(accelerate_bool),
               ## Epic-H WH-e: newton_kl TSVD truncation ratio (default 1e-8)
               as.double(newton_tsvd_ratio),
               ## Tikhonov ridge on dual λ (default 0.0 = off)
               as.double(ridge_lambda),
               PACKAGE = "leafblower")

  weights <- raw$weights
  calib_result    <- raw$result

  # WU-D: nest SOR diagnostics under $sor for clean namespace.
  # The C bridge always returns sor_min_omega and sor_n_damped as flat fields;
  # wrap them here so callers use result$sor$min_omega and result$sor$n_damped.
  calib_result$sor <- list(
    min_omega    = calib_result$sor_min_omega,
    n_damped     = calib_result$sor_n_damped,
    omega_mean   = calib_result$sor_omega_mean,
    any_latched  = calib_result$sor_any_latched,
    n_pinned_fb  = calib_result$sor_n_pinned_fb,
    n_warmup_fb  = calib_result$sor_n_warmup_fb,
    n_conv_fb    = calib_result$sor_n_conv_fb,
    n_resid_grew = calib_result$sor_n_resid_grew,
    n_monotone_cd = calib_result$sor_n_monotone_cd
  )
  calib_result$sor_min_omega   <- NULL
  calib_result$sor_n_damped    <- NULL
  calib_result$sor_omega_mean  <- NULL
  calib_result$sor_any_latched <- NULL
  calib_result$sor_n_pinned_fb <- NULL
  calib_result$sor_n_warmup_fb <- NULL
  calib_result$sor_n_conv_fb   <- NULL
  calib_result$sor_n_resid_grew  <- NULL
  calib_result$sor_n_monotone_cd <- NULL

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
      if (is.null(s) || is.na(s)) NA_character_
      else switch(as.character(s),
        "0" = "criterion",
        "4" = "budget",
        "3" = "error",
        "2" = "infeasible",
        "1" = "legacy",
        # leafblower-8eod: use solver-emitted convergence_stall_kind (set at RK_ERR_STALL
        # emission site) instead of the user-input accelerate_bool heuristic.
        # Audit (2026-05-08): raking.cpp:494 only fires inside the !st.accelerate branch
        # (raking.cpp:392 else-guard) → stall_kind=2. oris.cpp:188 fires for both SRAA
        # (accelerate=TRUE, stall_kind=1) and plain-BCD (accelerate=FALSE, stall_kind=2) —
        # NOT bijective with user flag; required route (a). stall_kind=0 → NA (no stall).
        "5" = {
          sk <- calib_result$convergence_stall_kind
          if (is.null(sk) || is.na(sk)) sk <- 0L
          switch(as.character(sk),
            "1" = "stall_wchange",
            "2" = "stall_kl",
            NA_character_
          )
        },
        NA_character_  # default
      )
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
  # NULL-guard required: `NULL == 2L` returns logical(0), and `if (logical(0))`
  # silently skips the branch — leaking garbage weights to caller.
  if (!is.null(calib_result$status) && calib_result$status == 2L)
    stop("leafblower: ", if (nchar(calib_result$message) > 0) calib_result$message
         else "infeasible problem")
  if (!is.null(calib_result$status) && calib_result$status == 3L)
    stop("leafblower: invalid arguments \u2014 ", calib_result$message)

  # Solver returns sum(weights) = n (enforced in src/oris.cpp, src/raking.cpp).
  # No wrapper-level
  # normalization — removing it preserves the bounds_mode="unit" strict-bounds
  # guarantee (oris's water-fill clamps are final; not re-pushed by post-scale).

  # NOTE: No post-normalization clamp to [min_weight, max_weight]. Clamping
  # here would break sum(weights * d) == target totals when per-cell mixing
  # parameters d are non-uniform: individual weights may legitimately exceed
  # per-cell bounds after expansion even when cell aggregates are in range.
  # The ORIS/raking solvers enforce bounds on the cell aggregate X[c], which
  # is the invariant that preserves calibration. See
  # tests/testthat/test-oris-nonuniform-d.R.

  if (calib_result$status == 4L) {
    e_final  <- calib_result$best_error
    b_iter   <- calib_result$best_iter
    iters    <- calib_result$iterations
    mstr     <- calib_result$convergence_used$metric
    if (is.na(mstr)) mstr <- "<metric>"
    tol_used <- if (conv$absolute_tol > 0) conv$absolute_tol else conv$pct_tol

    # Stall detection: if best was found much earlier than budget, solver is at a fixed point
    stall_ratio <- if (iters > 0L) b_iter / iters else 1.0
    if (is.numeric(stall_ratio) && stall_ratio < 0.5 &&
        is.numeric(e_final) && is.finite(e_final)) {
      warning(sprintf(
        "leafblower: fixed point at %s=%.2e (best at iter %d of %d, ratio=%.2f). More iterations will not improve calibration. Try: accelerate=TRUE, method='newton_kl', or method='oris+accel'.",
        mstr, e_final, b_iter, iters, stall_ratio),
        call. = FALSE)
    } else {
      # True budget: solver still converging when iterations ran out.
      # Use asymptotic rate from last check interval if available.
      e_prev    <- calib_result$metric_prev_check
      prev_iter <- calib_result$prev_check_iter
      interval  <- b_iter - prev_iter
      has_prev  <- is.numeric(e_prev) && is.finite(e_prev) &&
                   is.numeric(e_final) && is.finite(e_final) &&
                   e_prev > e_final && e_final > 0 &&
                   interval > 0L && is.numeric(tol_used) && tol_used > 0

      if (has_prev) {
        r_est  <- (e_final / e_prev)^(1 / interval)
        if (is.finite(r_est) && r_est > 0 && r_est < 1) {
          n_more  <- ceiling(log(tol_used / e_final) / log(r_est))
          n_total <- b_iter + n_more
          warning(sprintf(
            "leafblower: budget exhausted — %s=%.2e at %d iters. Asymptotic rate r=%.4f (last %d iters): ~%.0f total iterations needed.",
            mstr, e_final, iters, r_est, interval, n_total),
            call. = FALSE)
        } else {
          warning(sprintf(
            "leafblower: budget exhausted — %s=%.2e at %d iters. Increase max_iterations.",
            mstr, e_final, iters),
            call. = FALSE)
        }
      } else {
        warning(sprintf(
          "leafblower: budget exhausted — %s=%.2e at %d iters. Increase max_iterations.",
          mstr, e_final, iters),
          call. = FALSE)
      }
    }
  }
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
             "(max_weight=%.4g). Consider method='raking' or 'oris'."),
      calib_result$max_error, length(target), max_weight),
      call. = FALSE)
  }

  # Tier-1/2 calibration quality metrics.
  # Computed post-solver at O(nK) cost via helper: margin_kl (Tier 1),
  # weight_kl / design_effect / effective_observations (Tier 2).
  qm <- compute_quality_metrics(weights, target, data)
  calib_result[names(qm)] <- qm

  # --- c8w1: attach sparseness diagnostics to result ---
  calib_result$diagnostics <- list(
    sparseness = list(
      sparse_categories  = sparse_diag,
      pct_bounds_clamped = if (!is.null(calib_result$n_bounds_clamped))
                             calib_result$n_bounds_clamped / n_obs
                           else 0,
      thresholds         = list(target = 0.01, obs = 30L)
    )
  )

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

compute_sparseness_diag <- function(target, data, cat_threshold = 0.01, obs_threshold = 30L,
                                    na_margins = character(0)) {
  # Pre-hoist: extract data column names once; avoids repeated names(data) in loop.
  data_names  <- names(data)
  sparse_cats <- list()
  for (v in names(target)) {
    if (!v %in% data_names) next
    tv   <- target[[v]]           # hoist target[[v]] — one string lookup per margin
    col  <- data[[v]]             # hoist data[[v]]   — one string lookup per margin
    lvs  <- names(tv)             # hoist names() — reused for tabulate and inner loop
    # tabulate(match()) replaces table(): avoids factor()+unique() passes over col,
    # cutting 70% of per-call cost (profile: factor+unique.default = 66% self-time).
    if (v %in% na_margins) {
      # yaye: this margin has an injected "NA" bin. Mirror the solver encoding
      # (group_ids_r else-branch :475-477): conflate true-NA + literal-string-"NA"
      # into the "NA" level. match(NA, lvs) returns NA_integer_ and is dropped by
      # tabulate, so without the as.character()+fill the true-NA bin counts 0
      # (true-NA-blind) and is mis-flagged sparse on every add_na_proportion run.
      key  <- as.character(col)
      key[is.na(key)] <- "NA"
      cnts <- tabulate(match(key, lvs), nbins = length(lvs))
    } else {
      cnts <- tabulate(match(col, lvs), nbins = length(lvs))
    }
    names(cnts) <- lvs
    for (lv in lvs) {
      T_kj <- tv[[lv]]
      n_kj  <- cnts[[lv]]         # integer (tabulate output); 0 for absent levels
      if (T_kj < cat_threshold || n_kj < obs_threshold) {
        sparse_cats[[v]] <- c(sparse_cats[[v]], list(list(level = lv, T_kj = T_kj, n_kj = n_kj)))
      }
    }
  }
  sparse_cats
}

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

# CR-F8 (dtkn.8): exact-match validation for user-facing string arguments. base
# match.arg partial-matches (method="sink" -> "sinkhorn", metric="max" -> "max_err"),
# silently accepting typos/abbreviations, and emits an opaque "'arg' matches multiple
# formal arguments" error on ambiguous partials (method="or" -> oris/oris_soft).
# Require an exact choice, else raise a clear error naming the valid choices.
match_exact <- function(value, choices, arg_name) {
  if (length(value) != 1L || !is.character(value) || !(value %in% choices))
    stop(sprintf("%s must be exactly one of: %s; got: %s", arg_name,
                 paste0("\"", choices, "\"", collapse = ", "), deparse(value)),
         call. = FALSE)
  value
}

map_method <- function(method, verbose = 0) {
  method <- tolower(method)
  if (method %in% c("rake", "nrake")) {
    warning("method='", method, "' (IPF); using raking", call. = FALSE)
    method <- "raking"
  } else if (method == "nr") {
    warning("method='nr' (Newton-Raphson); using newton_kl", call. = FALSE)
    method <- "newton_kl"
  }
  match_exact(method, c("auto", "oris", "oris_soft", "raking", "sinkhorn",
                        "chebyshev", "greg", "greenkhorn", "logit", "newton_kl"),
              "method")
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
  # META.2: resolve pct alias to l1_weight before metric_int lookup
  if (!is.null(metric_raw) && identical(metric_raw, "pct")) metric_raw <- "l1_weight"
  metric <- match_exact(metric_raw,
    c("max_err", "mean_err", "kl", "chi2", "grake_norm", "l1_weight", "marginal_kl"),
    "convergence$metric")

  # pct is autumn/anesrake compatible: stops when Σ|Δw| STOPS IMPROVING (plateau).
  # Default its rule to "plateau" when pct is specified without an explicit rule.
  # absolute= with no explicit rule maps to "threshold" (hard stopping criterion).
  rule_explicit <- !is.null(convergence[["rule"]])
  rule_default  <- if (!rule_explicit && explicit_pct) "plateau"
                   else if (!rule_explicit && explicit_abs && !explicit_pct) "threshold"
                   else "improvement"
  rule_raw      <- convergence[["rule"]] %||% rule_default
  rule         <- match_exact(rule_raw, c("threshold", "improvement", "plateau"),
                              "convergence$rule")

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

  # CR-F9 (dtkn.9): the (0,1) plateau range check above guarded ONLY the explicit
  # `tol` entry path. The `pct` shorthand (convergence=list(pct=5)) resolves to the
  # same plateau pct_tol but bypassed it, silently admitting values the long form
  # rejects. Apply the identical range check to the resolved pct-shorthand tol.
  if (explicit_pct && rule == "plateau" && (pct_tol <= 0 || pct_tol >= 1))
    stop("convergence$pct must be in (0,1) for rule='plateau'")

  stop_when <- match_exact(convergence[["stop_when"]] %||% "any", c("any", "all"),
                           "convergence$stop_when")
  list(pct_tol = pct_tol, absolute_tol = absolute_tol,
       metric = metric, rule = rule, stop_when = stop_when)
}

parse_sor <- function(sor) {
  `%||%` <- function(a, b) if (is.null(a)) b else a
  if (is.null(sor)) {
    return(list(enabled = 0L, auto = 0L, omega_init = 1.0,
                omega_min = 0.3, omega_max = 1.5, omega_fixed = -1.0, burnin = 20L,
                omega_mode_id = 2L))
  }
  valid_keys <- c("auto", "omega_min", "omega_max", "omega", "omega_init", "burnin",
                  "omega_mode_id")
  bad <- setdiff(names(sor), valid_keys)
  if (length(bad))
    stop(sprintf("Unknown sor key(s): %s. Valid keys: %s",
                 paste(bad, collapse = ", "),
                 paste(valid_keys, collapse = ", ")))
  # omega_mode_id: 0=heuristic (0.7/1.05 nudge), 1=fixed (omega_max),
  #                2=iterate-change (free-coord ‖ΔX_free‖² estimator, default, e18t.9 SHIP).
  #                e18t.5 NO-SHIP rescinded: v2 iterate-change passes stepstone (50 vs 300 limit).
  raw_mode <- sor[["omega_mode_id"]]
  omega_mode_id <- if (is.null(raw_mode)) {
    2L
  } else if (is.character(raw_mode)) {
    switch(raw_mode,
      "heuristic" = 0L,
      "fixed"     = 1L,
      "spectral"  = 2L,
      stop(sprintf("Unknown omega_mode_id string '%s'. Use 'heuristic', 'fixed', or 'spectral'.",
                   raw_mode))
    )
  } else {
    as.integer(raw_mode)
  }
  list(
    enabled       = 1L,
    auto          = if (isTRUE(sor[["auto"]])) 1L else 0L,
    omega_init    = as.double(sor[["omega_init"]] %||% 1.0),
    omega_min     = as.double(sor[["omega_min"]] %||% 0.3),
    omega_max     = as.double(sor[["omega_max"]] %||% 1.5),
    omega_fixed   = as.double(sor[["omega"]] %||% -1.0),
    burnin        = as.integer(sor[["burnin"]] %||% 20L),
    omega_mode_id = omega_mode_id
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

compute_quality_metrics <- function(weights, target_list, df) {
  # Compute Tier-1/2 calibration quality metrics.
  # Returns list: design_effect, effective_observations, weight_kl, margin_kl.
  # Called post-solver at O(nK) cost from main harvest() body.
  if (!(length(weights) > 0L && length(target_list) > 0L &&
        is.finite(sum(weights)) && sum(weights) > 0)) {
    return(list(design_effect = NA_real_, effective_observations = NA_real_,
                weight_kl = NA_real_, margin_kl = NA_real_))
  }
  Z  <- sum(weights)
  n  <- length(weights)
  wn <- weights * n / Z                         # normalize to mean=1
  # leafblower-on7a v3 boundary: compute_quality_metrics uses Kish (1965) deff_K
  # (cheap observation-level path). H&V 2015 4-arg deff_H lives in design_effect(w, y, data, target).
  deff <- n * sum(weights^2) / Z^2

  # Per-margin KL helper: shared finalization across both single-pass and K-pass.
  # Inputs: T_k (named target probs), W_k (named observed probs aligned to a level set).
  margin_kl_one <- function(T_k, W_k) {
    if (any(T_k[setdiff(names(T_k), names(W_k))] > 0)) return(Inf)
    common <- intersect(names(T_k), names(W_k))
    T_sub  <- T_k[common]; W_sub <- W_k[common]
    pos    <- T_sub > 0                          # limit: 0 * log(0/W) = 0
    if (!any(pos)) return(0)
    sum(T_sub[pos] * log(T_sub[pos] / pmax(W_sub[pos], 1e-15)))
  }

  # B2 dispatch: single-pass cell-key aggregation when no NA in margin cols
  # AND K >= 3 (paste/rowsum overhead unprofitable for K < 3) AND all margin
  # cols are factors (precondition for integer-radix key encoding).
  margin_cols <- names(target_list)
  col_nlevels <- vapply(margin_cols, function(k) nlevels(df[[k]]), integer(1))
  use_single_pass <- length(target_list) >= 3L &&
    !anyNA(df[margin_cols]) &&
    all(vapply(margin_cols, function(k) is.factor(df[[k]]), logical(1)))

  # Guard: radix key must stay within 2^53 double-precision mantissa.
  # With K>=10 margins x ~50 levels, prod(nlevels) can exceed 2^53,
  # causing integer-key collisions and silent KL corruption.
  if (use_single_pass) {
    max_cells <- prod(col_nlevels)
    if (max_cells > 2^53) {
      use_single_pass <- FALSE
    }
  }

  margin_kl_value <- tryCatch({
    if (use_single_pass) {
      # Build integer cell-key by mixed-radix encoding of factor codes.
      # Faster than paste() for many-K, large-n: avoids string allocation.
      # rowsum() then aggregates weights to unique cells in one pass.
      key <- numeric(n)
      multiplier <- 1
      for (k in margin_cols) {
        key <- key + (as.integer(df[[k]]) - 1L) * multiplier
        multiplier <- multiplier * col_nlevels[[k]]
      }
      w_cell <- withCallingHandlers(
        rowsum(weights, key, na.action = NULL, reorder = TRUE),
        warning = function(w) {
          if (grepl("missing values", conditionMessage(w))) invokeRestart("muffleWarning")
        }
      )
      # Map each unique cell back to a representative row (for level lookup).
      cell_to_row <- match(as.numeric(rownames(w_cell)), key)
      vec_w_cell  <- as.vector(w_cell)
      sum(vapply(margin_cols, function(k) {
        T_k <- target_list[[k]]
        # Aggregate cell-level weights to margin-k levels.
        cell_levels_k <- df[[k]][cell_to_row]
        W_k <- tapply(vec_w_cell, cell_levels_k, sum) / Z
        margin_kl_one(T_k, W_k)
      }, numeric(1)))
    } else {
      # K-pass fallback: handles NA in margin cols (per-margin valid mask) and
      # K < 3 (where single-pass overhead is unprofitable).
      sum(sapply(margin_cols, function(k) {
        T_k   <- target_list[[k]]
        obs_k <- df[[k]]
        # CR-F1b (dtkn.12): under add_na_proportion=TRUE the target gains an "NA"
        # level (na_frac>0) and the solver encodes NA obs AS the "NA" level (the
        # NA->"NA" recode at ~line 507-516). Mirror that here: recode NA -> "NA"
        # and normalize over ALL obs, so W_k gains an "NA" mass matching the target
        # (else margin_kl_one sees a positive target "NA" with no W_k "NA" -> Inf).
        # SIGNAL = target level name "NA" present (all this helper receives). For the
        # in-scope add_na_proportion path this matches the solver encoding exactly.
        # LIMITATION (dtkn.13): a HAND-BUILT "NA" target level WITHOUT
        # add_na_proportion diverges — the solver factor path maps real NAs to gid
        # -1 (excluded), while this counts them into the "NA" bin. Hand-built "NA"
        # targets are ambiguous (harvest.R:141 doc); tracked separately, not fixed
        # here (this helper cannot see the injected-NA-margins set).
        if ("NA" %in% names(T_k)) {
          obs_chr <- as.character(obs_k)
          obs_chr[is.na(obs_k)] <- "NA"
          Z_k <- sum(weights)
          if (Z_k == 0) return(NA_real_)
          W_k <- tapply(weights, obs_chr, sum) / Z_k
          return(margin_kl_one(T_k, W_k))
        }
        valid <- !is.na(obs_k)                  # exclude NA-coded observations
        w_v   <- weights[valid]
        Z_k   <- sum(w_v)
        if (Z_k == 0) return(NA_real_)
        obs_kv <- obs_k[valid]
        # droplevels() has no character method — guard to factor; tapply coerces char/numeric/logical itself [mb06]
        grp    <- if (is.factor(obs_kv)) droplevels(obs_kv) else obs_kv
        W_k    <- tapply(w_v, grp, sum) / Z_k
        margin_kl_one(T_k, W_k)
      }), na.rm = FALSE)
    }
  }, error = function(e) { warning("margin_kl: ", conditionMessage(e), call. = FALSE); NA_real_ })

  list(
    design_effect = deff,
    effective_observations = if (is.finite(deff) && deff > 0) n / deff else NA_real_,
    weight_kl = sum(wn * log(pmax(wn, 1e-15))) / n,
    margin_kl = margin_kl_value
  )
}
