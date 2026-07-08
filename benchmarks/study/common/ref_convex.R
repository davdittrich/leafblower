# benchmarks/study/common/ref_convex.R
#
# Independent high-precision convex reference solver -- the RQ5 agreement
# anchor (DESIGN.md Section 6 / Blocker G). Ticket leafblower-2ouc.5 (WU-4).
#
# INDEPENDENCE (Mechanism/Forbidden, DESIGN.md): this file calls CVXR only.
# It NEVER calls leafblower::harvest() or any leafblower method -- it is not
# a benchmarked competitor, it is the correctness anchor other solvers (incl.
# leafblower's own) are measured against. The per-family objective below is
# the textbook distance function for that family (not an arbitrary
# alternative) so the anchor targets the SAME mathematical problem
# leafblower's own solvers target; only the numerical solve is independent.
#
# Family objectives (derivations cited against leafblower's OWN solver
# source, read-only, for definition purposes only -- never invoked):
#   kl     Sum w_i*log(w_i/d_i)                       (src/raking.cpp: cyclic
#          KL projections converge to "the bounded KL minimum", Csiszar-
#          Tusnady 1984). CVXR: sum(kl_div(w,d)) == this objective + a
#          constant on the feasible set (kl_div(x,y)=x*log(x/y)-x+y; the
#          margin constraints below pin Sum(w)==Sum(d)==N so the -w+d term
#          is constant at every feasible point).
#   chi2   Sum (w_i-d_i)^2/d_i                        (Deville-Sarndal
#          "linear"/chi-square GREG distance; DESIGN.md Section 2 "chi2
#          linear -- GREG").
#   logit  Sum[(w_i-L)log(w_i-L) + (U-w_i)log(U-w_i)] (re-derived from
#          leafblower's ACTUAL logit KKT structure: src/logit_calib.cpp:286-
#          296 solves w_c = L_c + (U_c-L_c)*sigma(z_c) with z_c a linear
#          combination of duals -- i.e. dD/dw = log((w-L)/(U-w)) = z at the
#          optimum. Integrating gives the objective above (verified:
#          d/dw[(w-L)log(w-L)+(U-w)log(U-w)] = log((w-L)/(U-w))). NOTE this
#          is a symmetric binary-entropy-type distance with NO d_i anchor
#          term -- logit_calib.cpp's design-weight use is confined to a
#          Newton warm-start (Layer 2), not the objective/halt test. This
#          differs from the classical Deville-Sarndal(1992) asymmetric-scale
#          logit link in src/logit.hpp (LinkFn, unused by logit_calib.cpp);
#          the anchor here matches the solver actually shipped.
#   minimax  min_w max_{k,j} |S_kj/N - T_kj|  s.t. Sum(w)==N (population
#          total preserved -- standard calibration convention; also the
#          convention that makes this genuinely an LP, matching DESIGN.md
#          Section 2's "Py: cvxpy LP baseline" framing, rather than the
#          linear-fractional S/W leafblower's chebyshev.cpp IPM tracks with
#          a floating W). Blocker G: the L-inf optimum is a non-unique face,
#          so only the ACHIEVED OBJECTIVE VALUE is a valid anchor, never the
#          weight vector -- solve_ref() never returns mode="weight_vector"
#          for minimax, and the objective value stored is recomputed via
#          metrics.R::margin_stats() (the SAME golden-checked code path
#          every other solver's rows are scored with), not the raw LP `t`.
#
# Scope limit (DESIGN.md Section 6, ticket Section VI): a 1e-12 convex solve
# is infeasible at stepstone scale (n ~ 1.58M) -- REF_MAX_N below refuses
# any problem larger than that, loudly, rather than silently hanging or
# faking a result. 1.58M has NO anchor; this is a stated limitation, not a
# bug to work around.

suppressPackageStartupMessages({
  library(CVXR)
})

.rc_script_dir <- function() {
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  if (length(f)) dirname(f) else "benchmarks/study/common"
}
source(file.path(.rc_script_dir(), "metrics.R"))

REF_SOLVER_TOL <- 1e-12
REF_MAX_N <- 50000L
REF_STRICTLY_CONVEX_FAMILIES <- c("kl", "chi2", "logit")
REF_FAMILIES <- c("kl", "chi2", "logit", "minimax")
REF_ANCHOR_KIND <- list(kl = "weight_vector", chi2 = "weight_vector",
                         logit = "weight_vector", minimax = "objective_value")

.rc_assert_scope <- function(problem) {
  n <- nrow(problem$data)
  if (n > REF_MAX_N) {
    stop("ref_convex: problem '", problem$id, "' has n=", n, " > REF_MAX_N=",
         REF_MAX_N, " -- a ~1e-12 convex reference solve is infeasible at ",
         "this scale (DESIGN.md Section 6: '1.58M has no anchor', stated ",
         "limitation, never faked). Exclude from the RQ5 anchor.", call. = FALSE)
  }
}

# Population total N and per-margin category target COUNTS (T_kj * N), plus
# the per-observation category label vector for each margin. Shared by every
# family: the margin constraint is always Sum_{i in cat kj} w_i == T_kj * N,
# a LINEAR constraint on w regardless of family (only the objective differs).
.rc_margin_targets <- function(problem) {
  N <- sum(problem$design_weights)
  setNames(lapply(problem$margins, function(m) {
    list(groups = as.character(problem$data[[m]]),
         target_counts = problem$targets[[m]] * N)
  }), problem$margins)
}

.rc_objective <- function(family, w, d, L, U) {
  if (family == "kl") {
    if (any(d <= 0)) {
      stop("ref_convex: kl family requires d_i > 0 for all i (got d_i<=0); ",
           "exclude zero-design-weight rows before calling solve_ref().", call. = FALSE)
    }
    sum(kl_div(w, d))
  } else if (family == "chi2") {
    if (any(d <= 0)) {
      stop("ref_convex: chi2 family requires d_i > 0 for all i (got d_i<=0); ",
           "exclude zero-design-weight rows before calling solve_ref().", call. = FALSE)
    }
    sum(square(w - d) * (1 / d))
  } else if (family == "logit") {
    if (!is.finite(U)) {
      stop("ref_convex: logit family requires a finite max bound (matches ",
           "leafblower's own guard, src/logit_calib.cpp:46 'max_weight must ",
           "be finite and positive').", call. = FALSE)
    }
    -sum(entr(w - L)) - sum(entr(U - w))
  } else {
    stop("ref_convex: unsupported strictly-convex family: ", family, call. = FALSE)
  }
}

#' Solve the independent convex reference problem for one (problem, family).
#'
#' Strictly-convex families (kl, chi2, logit) return mode="weight_vector"
#' (unique optimum). minimax returns mode="objective_value" ONLY (Blocker G:
#' non-unique optimum face) -- weights are never exposed as a stored anchor
#' for minimax, only used internally to recompute the achieved L-inf error
#' via the shared golden metrics.R::margin_stats().
#'
#' @param problem standardized problem object (problem_io.R::load_problem_spec())
#' @param family one of REF_FAMILIES
#' @param tol CLARABEL gap/feasibility tolerance (default REF_SOLVER_TOL)
solve_ref <- function(problem, family, tol = REF_SOLVER_TOL) {
  if (!family %in% REF_FAMILIES) {
    stop("ref_convex: unknown family '", family, "' (expected one of: ",
         paste(REF_FAMILIES, collapse = ", "), ")", call. = FALSE)
  }
  .rc_assert_scope(problem)

  n <- nrow(problem$data)
  d <- problem$design_weights
  L <- problem$bounds$min
  U <- problem$bounds$max
  N <- sum(d)
  mt <- .rc_margin_targets(problem)

  w <- Variable(n)
  base_cons <- list(w >= L, sum(w) == N)
  if (is.finite(U)) base_cons <- c(base_cons, list(w <= U))

  # Built via explicit accumulation (not nested lapply()+unlist(recursive=FALSE))
  # because unlist() on a list of CVXR S4 Constraint objects can silently
  # mis-flatten nested list-of-list structures (see minimax linf_cons below,
  # which hit exactly this bug with 2 constraints per level).
  eq_cons <- list()
  for (mm in mt) {
    lv <- names(mm$target_counts)
    for (l in lv) {
      eq_cons[[length(eq_cons) + 1L]] <- sum(w[mm$groups == l]) == unname(mm$target_counts[[l]])
    }
  }

  solve_ctl <- list(solver = "CLARABEL", tol_gap_abs = tol, tol_gap_rel = tol,
                     tol_feas = tol, tol_infeas_abs = tol, tol_infeas_rel = tol)

  if (family %in% REF_STRICTLY_CONVEX_FAMILIES) {
    obj <- .rc_objective(family, w, d, L, U)
    prob <- Problem(Minimize(obj), c(base_cons, eq_cons))
    do.call(psolve, c(list(problem = prob), solve_ctl))
    st <- status(prob)
    if (!st %in% c("optimal", "optimal_inaccurate")) {
      stop("ref_convex: CLARABEL did not reach optimality (family=", family,
           ", problem=", problem$id, "): status=", st, call. = FALSE)
    }
    list(family = family, mode = "weight_vector", weights = as.numeric(value(w)),
         objective = value(prob), solver_status = st)
  } else { # minimax
    t <- Variable(1)
    linf_cons <- list()
    for (mm in mt) {
      lv <- names(mm$target_counts)
      for (l in lv) {
        S <- sum(w[mm$groups == l]); Tc <- unname(mm$target_counts[[l]])
        linf_cons[[length(linf_cons) + 1L]] <- S - Tc <= t * N
        linf_cons[[length(linf_cons) + 1L]] <- Tc - S <= t * N
      }
    }
    prob <- Problem(Minimize(t), c(base_cons, list(t >= 0), linf_cons))
    do.call(psolve, c(list(problem = prob), solve_ctl))
    st <- status(prob)
    if (!st %in% c("optimal", "optimal_inaccurate")) {
      stop("ref_convex: CLARABEL did not reach optimality (family=minimax, ",
           "problem=", problem$id, "): status=", st, call. = FALSE)
    }
    w_solved <- as.numeric(value(w))
    groups <- setNames(lapply(problem$margins, function(m) as.character(problem$data[[m]])),
                        problem$margins)
    ms <- margin_stats(w_solved, groups, problem$targets)
    list(family = "minimax", mode = "objective_value", obj_val = ms$margin_linf,
         weights = w_solved, solver_status = st)
  }
}

ref_weights_path <- function(family, problem_id, out_dir = "benchmarks/study/results") {
  file.path(out_dir, "weights", sprintf("ref_%s__%s__t1__na.parquet", family, problem_id))
}

ref_objective_path <- function(family, problem_id, out_dir = "benchmarks/study/results") {
  file.path(out_dir, "ref_objective", sprintf("%s__%s.json", family, problem_id))
}

#' Store a solve_ref() result as a pseudo-solver row (strictly-convex
#' families: weights/ parquet, one column `weight`, length n) or the
#' achieved-objective-value record (minimax: small JSON, per contract.md
#' weights_ref-adjacent convention -- never a weight-vector parquet, so
#' minimax can never be accidentally consumed by RQ5 vector correlation).
store_ref <- function(problem, family, result, out_dir = "benchmarks/study/results") {
  if (result$mode == "weight_vector") {
    path <- ref_weights_path(family, problem$id, out_dir)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    arrow::write_parquet(data.frame(weight = result$weights), path)
    path
  } else {
    path <- ref_objective_path(family, problem$id, out_dir)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    jsonlite::write_json(list(problem = problem$id, family = family,
                               anchor_kind = "objective_value",
                               obj_val = result$obj_val,
                               solver_status = result$solver_status,
                               note = "Blocker G: L-inf optimum is a non-unique face; excluded from weight-vector correlation (DESIGN.md Section 6)."),
                          path, auto_unbox = TRUE, digits = 15)
    path
  }
}
