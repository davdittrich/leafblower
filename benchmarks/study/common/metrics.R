# benchmarks/study/common/metrics.R
#
# Shared quality metrics for the leafblower benchmark study (DESIGN.md §6).
# Reconciles THREE existing, mutually-inconsistent implementations:
#   benchmarks/stepstone_all_methods.R:22   fit_metrics      (:34 t>0&S>0 gate BUG, :37 d_i=1 BUG)
#   benchmarks/python_ipf_benchmark.py:31   compute_metrics  (same two bugs)
#   benchmarks/allmethod_bench.R:30         compute_metrics  (3rd, differently-normalised weight_kl/chi2)
# This file is the canonical, bug-fixed replacement. Read-only w.r.t. the three
# files above — none of them are modified here.
#
# Guards enforced (DESIGN.md §6, ticket leafblower-2ouc.2):
#  - margin KL is DIVERGENT-AWARE: a starved category (T_kj>0, p_kj->0) yields
#    +Inf and sets a `divergent` flag. The term is NEVER dropped (no T>0&p>0
#    gate) — that gate is precisely the bug that flatters catastrophic solver
#    failures in the two `t>0 & S>0`-gated legacy impls above.
#  - weight_kl uses the REAL per-problem design weights d_i (never hardcoded
#    d_i=1) and is reported family-native with an explicit neutral-axis
#    caveat: it IS the KL/raking dual objective, so it flatters that family.
#  - DEFF is reported as "Kish weighting DEFF/UWE" (raw-w) — NOT the true
#    design effect — plus a SEPARATE g-weight efficiency 1+CV^2(g), g=w/d,
#    whenever d_i != 1 (Deville-Saerndal; conflated base-design and
#    calibration-injected variance otherwise).
#  - No cancellation: DEFF/ESS/CV^2 use ratio-of-sums-of-squares forms
#    (n*Sum(x^2)/Sum(x)^2), never Mean(x^2)-Mean(x)^2. Zero-design-weight
#    rows (d_i==0) are excluded from weight_kl and g-weight stats (0*Inf
#    guard on a divergent d_i=0 shift).
#  - RQ5 agreement: weight-vector Pearson/Spearman/max|dw|/cosine ONLY for
#    strictly-convex families (kl, chi2, logit); minimax is judged on
#    achieved-L-infinity OBJECTIVE VALUE agreement (Blocker G — the L-inf LP
#    optimum lies on a face, not a unique vertex, so weight vectors from two
#    equally-correct minimax solves legitimately differ).

`%||%` <- function(a, b) if (!is.null(a)) a else b

KL_STRICTLY_CONVEX_FAMILIES <- c("kl", "chi2", "logit")
KL_NATIVE_FAMILIES <- c("kl", "raking", "sinkhorn", "greenkhorn", "oris", "oris_soft", "newton_kl")

# ---- margin KL / L-inf / L1 -------------------------------------------------

#' One margin's proportions, KL term, and error stats.
#' @param w numeric weight vector, length n
#' @param group factor/character vector, length n, category per obs for this margin
#' @param target named numeric vector of target proportions, names = levels (sum to 1)
margin_one <- function(w, group, target) {
  W  <- sum(w)
  lv <- names(target)
  S  <- vapply(lv, function(l) sum(w[group == l]), numeric(1))
  p  <- unname(S / W)
  T  <- unname(target[lv])
  err <- abs(p - T)
  # Divergent-aware KL term: T*log(T/p). T==0 contributes 0 (no target mass to
  # miss). T>0 & p<=0 (starved category) -> +Inf, NEVER dropped.
  term <- ifelse(T == 0, 0,
            ifelse(p <= 0, Inf, T * log(T / p)))
  list(kl = sum(term), divergent = any(T > 0 & p <= 0),
       linf = max(err), l1 = sum(err), p = p, T = T)
}

#' Margin KL (mean/max across margins), divergence flag, margin L-inf/L1.
#' @param w numeric weight vector
#' @param groups named list of factor/character vectors, one per margin
#' @param targets named list of target-proportion vectors, same names as groups
margin_stats <- function(w, groups, targets) {
  ks  <- names(targets)
  per <- lapply(ks, function(k) margin_one(w, groups[[k]], targets[[k]]))
  names(per) <- ks
  kl_k      <- vapply(per, `[[`, numeric(1), "kl")
  divergent <- vapply(per, `[[`, logical(1), "divergent")
  linf_k    <- vapply(per, `[[`, numeric(1), "linf")
  l1_k      <- vapply(per, `[[`, numeric(1), "l1")
  list(
    marg_kl_mean = mean(kl_k), marg_kl_max = max(kl_k),
    divergent = any(divergent), divergent_margins = ks[divergent],
    margin_linf = max(linf_k), margin_l1 = sum(l1_k),
    per_margin = per
  )
}

# ---- ESS / Kish DEFF / g-weight efficiency ---------------------------------

# Ratio-of-sums-of-squares form of 1+CV^2(x) == n*Sum(x^2)/Sum(x)^2.
# Cancellation-free: no Mean(x^2)-Mean(x)^2 subtraction of close numbers.
.deff_ratio <- function(x) length(x) * sum(x^2) / sum(x)^2

#' ESS (Kish), Kish weighting DEFF/UWE (raw w), and g-weight efficiency
#' 1+CV^2(g) (g=w/d) when design weights are non-uniform.
#' @param w numeric weight vector
#' @param d design weights (d_i), same length as w; NULL => no g-weight stats
ess_deff <- function(w, d = NULL) {
  n   <- length(w)
  ess <- sum(w)^2 / sum(w^2)
  out <- list(
    ESS = ess, n = n,
    DEFF_kish = n / ess,   # == .deff_ratio(w); n/ESS is the design-quoted form
    DEFF_kish_label = "Kish weighting DEFF/UWE (raw w) -- NOT the true design effect",
    g_weighted = FALSE, DEFF_g = NA_real_
  )
  if (!is.null(d)) {
    use <- d > 0
    if (any(!use)) out$excluded_zero_design <- sum(!use)
    if (any(d[use] != 1)) {
      g <- w[use] / d[use]
      out$DEFF_g <- .deff_ratio(g)
      out$g_weighted <- TRUE
    }
  }
  out
}

# ---- closeness-to-design weight_kl (family-native) -------------------------

#' weight_kl = Sum w_i*log(w_i/d_i), REAL per-problem d_i (never hardcoded 1).
#' Rows with d_i==0 excluded (0*Inf guard on a divergent shift); w_i==0
#' contributes 0 by convention (0*log(0/d)=0).
#' @param family objective family string ("kl","chi2","logit","minimax",...);
#'   used only to attach the neutral-axis caveat, not the computation.
weight_kl <- function(w, d, family = NA_character_) {
  stopifnot(length(w) == length(d))
  use <- d > 0
  w_u <- w[use]; d_u <- d[use]
  term <- ifelse(w_u > 0, w_u * log(w_u / d_u), 0)
  native <- !is.na(family) && family %in% KL_NATIVE_FAMILIES
  note <- if (native)
    "weight_kl IS this family's objective (KL/raking dual) -- not a neutral cross-family axis"
  else
    "weight_kl reported as closeness-to-design diagnostic (non-native for this family)"
  list(weight_kl = sum(term), excluded_zero_design = sum(!use),
       family = family, neutral_axis_note = note)
}

# ---- bound violation --------------------------------------------------------

#' Count and max/mean magnitude of max(0, L-w, w-U). L/U scalar or per-obs vector.
bound_violation <- function(w, L = 0, U = Inf) {
  viol <- pmax(0, L - w, w - U)
  list(count = sum(viol > 0), max = max(viol), mean = mean(viol))
}

# ---- RQ5 agreement (Blocker G: minimax excluded from vector correlation) --

#' Weight-vector agreement for strictly-convex families; achieved-L-inf
#' objective-value agreement for minimax (unique-optimum caveat, Blocker G).
#' @param w,w_ref weight vectors under comparison (candidate, anchor)
#' @param family objective family string
#' @param obj_val,obj_val_ref achieved L-inf objective values (minimax only)
#' @param obj_tol absolute tolerance for objective-value agreement
agreement <- function(w, w_ref, family, obj_val = NULL, obj_val_ref = NULL, obj_tol = 1e-6) {
  if (!is.na(family) && family %in% KL_STRICTLY_CONVEX_FAMILIES) {
    list(family = family, mode = "weight_vector",
         pearson = stats::cor(w, w_ref, method = "pearson"),
         spearman = stats::cor(w, w_ref, method = "spearman"),
         max_abs_diff = max(abs(w - w_ref)),
         cosine = sum(w * w_ref) / (sqrt(sum(w^2)) * sqrt(sum(w_ref^2))))
  } else {
    stopifnot(!is.null(obj_val), !is.null(obj_val_ref))
    d <- abs(obj_val - obj_val_ref)
    list(family = family, mode = "objective_value",
         obj_val = obj_val, obj_val_ref = obj_val_ref,
         abs_diff = d, agree = d <= obj_tol)
  }
}

# ---- top-level convenience wrapper -----------------------------------------

#' Full metrics record for one (weights, problem) pair.
#' @param w numeric weight vector
#' @param groups named list of per-margin category vectors
#' @param targets named list of per-margin target-proportion vectors
#' @param d design weights (default: rep(1, n) -> d_i=1 uniform-start problems)
#' @param bounds list(L=, U=) or NULL
#' @param family objective family string (see KL_STRICTLY_CONVEX_FAMILIES / KL_NATIVE_FAMILIES)
#' @param w_ref,obj_val,obj_val_ref,obj_tol RQ5 agreement inputs (all optional)
compute_metrics <- function(w, groups, targets, d = NULL, bounds = NULL,
                             family = NA_character_,
                             w_ref = NULL, obj_val = NULL, obj_val_ref = NULL,
                             obj_tol = 1e-6) {
  if (is.null(d)) d <- rep(1, length(w))
  ms <- margin_stats(w, groups, targets)
  ed <- ess_deff(w, d)
  wk <- weight_kl(w, d, family)
  out <- list(
    n = length(w), W = sum(w),
    marg_kl_mean = ms$marg_kl_mean, marg_kl_max = ms$marg_kl_max,
    marg_kl_divergent = ms$divergent, marg_kl_divergent_margins = ms$divergent_margins,
    margin_linf = ms$margin_linf, margin_l1 = ms$margin_l1,
    ESS = ed$ESS, DEFF_kish = ed$DEFF_kish, DEFF_kish_label = ed$DEFF_kish_label,
    DEFF_g = ed$DEFF_g, g_weighted = ed$g_weighted,
    weight_kl = wk$weight_kl, weight_kl_excluded_zero_design = wk$excluded_zero_design,
    weight_kl_note = wk$neutral_axis_note,
    wmin = min(w), wmed = stats::median(w), wmax = max(w)
  )
  if (!is.null(bounds)) {
    bv <- bound_violation(w, bounds$L %||% 0, bounds$U %||% Inf)
    out$bound_viol_count <- bv$count; out$bound_viol_max <- bv$max; out$bound_viol_mean <- bv$mean
  }
  if (!is.null(w_ref)) {
    out$agreement <- agreement(w, w_ref, family, obj_val, obj_val_ref, obj_tol)
  }
  out
}
