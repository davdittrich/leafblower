# STUDY-BRANCH-ONLY-DO-NOT-MERGE
# replicate_weights.R -- Rao-Wu-Yue (RWY) rescaled bootstrap replicate
# design-weights for variance estimation over synthetic i.i.d. instances
# (WU leafblower-2ouc.14.1). Obs-as-PSU: no strata/clusters, since study
# instances are generated i.i.d. draws, not a stratified/clustered design.
#
# Reference: Rao, J.N.K. and Wu, C.F.J. (1988), "Resampling Inference with
# Complex Survey Data", JASA 83(401):231-241 (rescaled bootstrap).

#' Build Rao-Wu rescaled-bootstrap replicate design-weights.
#'
#' For each replicate b = 1..B, draws resample counts
#' r ~ Multinomial(m; 1/n, ..., 1/n) with m = n - 1, then forms the
#' mean-preserving rescale factor
#'   a_i = 1 - sqrt(m/(n-1)) + sqrt(m/(n-1)) * (n/m) * r_i
#' which, at m = n - 1, collapses to a_i = (n/(n-1)) * r_i (the first two
#' terms cancel exactly). Replicate weights are w_i^(b) = design_weights_i * a_i.
#' E[r_i] = m/n = (n-1)/n, so E[a_i] = (n/(n-1)) * (n-1)/n = 1: replicate
#' weights are mean-preserving in expectation. Some a_i may legitimately be
#' zero (unit not resampled in that replicate) -- this is standard bootstrap
#' behaviour and is NOT floored to a positive value.
#'
#' @param design_weights Numeric vector of per-observation design weights,
#'   length n >= 2, all finite and > 0.
#' @param B Integer number of bootstrap replicates, B >= 1.
#' @param seed Integer RNG seed. Same seed => byte-identical output matrix.
#'   The caller's global RNG state (.Random.seed) is saved and restored, so
#'   this call does not perturb the caller's RNG stream.
#' @param type Replicate-weight scheme. Only "RWY" (Rao-Wu-Yue rescaled
#'   bootstrap) is implemented; any other value raises an error.
#' @return An n x B numeric matrix of replicate design-weights.
rm_make_replicate_weights <- function(design_weights, B, seed, type = "RWY") {
  if (!identical(type, "RWY")) {
    stop(sprintf("rm_make_replicate_weights: unsupported type '%s' (only \"RWY\" is implemented)", type))
  }
  if (!is.numeric(design_weights)) {
    stop("rm_make_replicate_weights: design_weights must be numeric")
  }
  n <- length(design_weights)
  if (n < 2L) {
    stop("rm_make_replicate_weights: design_weights must have length n >= 2")
  }
  if (!all(is.finite(design_weights))) {
    stop("rm_make_replicate_weights: design_weights must all be finite")
  }
  if (!all(design_weights > 0)) {
    stop("rm_make_replicate_weights: design_weights must all be strictly positive")
  }
  if (!is.numeric(B) || length(B) != 1L || is.na(B) || B < 1 || B != as.integer(B)) {
    stop("rm_make_replicate_weights: B must be a single integer >= 1")
  }
  B <- as.integer(B)

  m <- n - 1L
  scale_a <- n / m  # = sqrt(m/(n-1)) * (n/m) since m == n-1 exactly

  draw <- function() {
    r <- stats::rmultinom(1L, size = m, prob = rep(1 / n, n))
    a <- scale_a * r[, 1L]
    design_weights * a
  }

  run <- function() {
    cols <- replicate(B, draw(), simplify = FALSE)
    matrix(unlist(cols), nrow = n, ncol = B)
  }

  if (requireNamespace("withr", quietly = TRUE)) {
    withr::with_seed(seed, run())
  } else {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
      } else {
        assign(".Random.seed", old_seed, envir = .GlobalEnv)
      }
    }, add = TRUE)
    set.seed(seed)
    run()
  }
}
