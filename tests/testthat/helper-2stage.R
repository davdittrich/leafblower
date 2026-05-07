# helper-2stage.R — DGP generators for 2-stage raking tests.
# In-test seeded helpers; no .rds files persisted.

#' Generate K=9 binary survey data from a Dirichlet(0.1) sparse joint.
#'
#' Amended DGP (spec §8, iter 3): N=10000; K=9 binary margins; Dirichlet(α=0.1)
#' joint; margin TARGETS skewed at target_skew (default 0.95) for category "1"
#' (i.e., 0.05/0.95 split), i.i.d. across all K margins; coarse = first 3
#' margins; min_cell_n = 30.
#'
#' With skewed targets (0.95) combined with max_weight=2 (tight bound regime),
#' single-stage raking must upweight the rare "0" observations by factors far
#' exceeding 2, causing bound violations on ≥80% of seeds. The 2-stage path
#' calibrates within coarse cells where local marginals are less extreme,
#' staying within bounds and achieving Σ|w·X_k/N − target_k| ≤ 1e-4.
#'
#' @param seed        Integer seed.
#' @param N           Sample size (default 10000, per spec §8).
#' @param K           Number of binary margins (default 9).
#' @param target_skew Probability of category "1" in targets (default 0.95).
#' @return list(df, targets, coarse_mask)
make_k9_sparse <- function(seed, N = 10000L, K = 9L, target_skew = 0.95) {
  set.seed(seed)
  n_cells <- 2L ^ K
  alpha   <- rep(0.1, n_cells)
  probs   <- rgamma(n_cells, shape = alpha)
  probs   <- probs / sum(probs)
  cells   <- sample(seq_len(n_cells), N, replace = TRUE, prob = probs)
  # Decode cell index to K binary columns (0/1).
  mat <- matrix(
    as.integer((outer(cells - 1L, 2L ^ (0L:(K - 1L)), bitwAnd) > 0L)),
    nrow = N, ncol = K,
    dimnames = list(NULL, paste0("x", seq_len(K)))
  )
  df <- as.data.frame(mat)
  for (j in seq_len(K)) df[[j]] <- as.factor(df[[j]])
  # Skewed targets: each margin targets target_skew for category "1".
  targets <- lapply(seq_len(K), function(k) {
    c(`0` = 1.0 - target_skew, `1` = target_skew)
  })
  names(targets) <- names(df)
  # Coarse = first 3 margins.
  coarse_mask <- c(rep(1L, 3L), rep(0L, K - 3L))
  list(df = df, targets = targets, coarse_mask = coarse_mask)
}

#' Generate rescue DGP for spec §8 v3 (K=6, N=80, chain-correlated binary).
#'
#' Spec §8 v3: N=80; K=6 binary margins (3 coarse + 3 fine); joint as chain
#' of skewed correlated margins. coarse_1 ~ Bern(0.1), coarse_2 ~ Bern(0.15),
#' coarse_3 ~ Bern(0.2); fine_k = ifelse(coarse_k==0, Bern(0.05), Bern(0.95)).
#' Targets = uniform c('0'=0.5, '1'=0.5) for every margin.
#' coarse_mask = c(1,1,1,0,0,0); min_cell_n = 30; default max_weight.
#'
#' @param seed Integer seed.
#' @param N    Sample size (default 80, per spec §8 v3).
#' @return list(df, targets, coarse_mask)
make_rescue_dgp <- function(seed, N = 80L) {
  set.seed(seed)
  # Coarse margins: skewed Bernoulli as specified.
  g1 <- rbinom(N, 1L, 0.10)
  g2 <- rbinom(N, 1L, 0.15)
  g3 <- rbinom(N, 1L, 0.20)
  # Fine margins: chain-correlated with their coarse counterpart.
  f1 <- ifelse(g1 == 0L, rbinom(N, 1L, 0.05), rbinom(N, 1L, 0.95))
  f2 <- ifelse(g2 == 0L, rbinom(N, 1L, 0.05), rbinom(N, 1L, 0.95))
  f3 <- ifelse(g3 == 0L, rbinom(N, 1L, 0.05), rbinom(N, 1L, 0.95))
  df <- data.frame(
    g1 = factor(g1), g2 = factor(g2), g3 = factor(g3),
    f1 = factor(f1), f2 = factor(f2), f3 = factor(f3)
  )
  # Uniform targets: 50/50 for every margin.
  targets <- lapply(names(df), function(nm) c(`0` = 0.5, `1` = 0.5))
  names(targets) <- names(df)
  coarse_mask <- c(1L, 1L, 1L, 0L, 0L, 0L)
  list(df = df, targets = targets, coarse_mask = coarse_mask)
}

#' Generate a weight-stress DGP that causes single-stage raking to fail.
#'
#' Uses extreme skew (10:1 ratio) in the first margin combined with
#' highly non-uniform fine margins, small N, and tight targets, such that
#' single-stage raking cannot satisfy all margins simultaneously and
#' diverges or returns NaN on most seeds.
#'
#' @param seed Integer seed.
#' @param N    Sample size (default 80).
#' @return list(df, targets, coarse_mask)
make_stress_dgp <- function(seed, N = 80L) {
  set.seed(seed)
  K <- 6L
  # Skewed coarse margins.
  grp1 <- sample(c(0L, 1L), N, replace = TRUE, prob = c(0.9, 0.1))
  grp2 <- sample(c(0L, 1L), N, replace = TRUE, prob = c(0.85, 0.15))
  grp3 <- sample(c(0L, 1L), N, replace = TRUE, prob = c(0.8, 0.2))
  # Fine margins with cross-correlations.
  fine1 <- ifelse(grp1 == 0L, sample(0:1, N, TRUE, c(0.95, 0.05)),
                              sample(0:1, N, TRUE, c(0.05, 0.95)))
  fine2 <- ifelse(grp2 == 0L, sample(0:1, N, TRUE, c(0.95, 0.05)),
                              sample(0:1, N, TRUE, c(0.05, 0.95)))
  fine3 <- ifelse(grp3 == 0L, sample(0:1, N, TRUE, c(0.9, 0.1)),
                              sample(0:1, N, TRUE, c(0.1, 0.9)))
  df <- data.frame(
    g1 = factor(grp1), g2 = factor(grp2), g3 = factor(grp3),
    f1 = factor(fine1), f2 = factor(fine2), f3 = factor(fine3)
  )
  # Targets that conflict with the skewed joint: require equal splits.
  targets <- lapply(names(df), function(nm) c(`0` = 0.5, `1` = 0.5))
  names(targets) <- names(df)
  coarse_mask <- c(1L, 1L, 1L, 0L, 0L, 0L)
  list(df = df, targets = targets, coarse_mask = coarse_mask)
}

#' Check margin residuals for a calibrated result.
#' For binary (0/1) margins: checks target proportion for category "1".
#' Assumes all margin targets are c(`0`=p0, `1`=p1) style.
#' For the rescue DGP, all targets are 0.5/0.5 so tgt_list defaults to NULL
#' and the check uses 0.5.
#' @return max |sum(w * X_k) / N - target_1_k| over all K margins.
max_margin_resid <- function(w, df, tgt_list = NULL) {
  K  <- ncol(df)
  N  <- length(w)
  mx <- 0
  for (k in seq_len(K)) {
    x     <- as.numeric(as.character(df[[k]]))
    tgt_1 <- if (!is.null(tgt_list)) {
      nm <- names(df)[k]
      unname(tgt_list[[nm]]["1"])
    } else {
      0.5
    }
    if (is.na(tgt_1)) tgt_1 <- 0.5
    err <- abs(sum(w * x) / N - tgt_1)
    if (err > mx) mx <- err
  }
  mx
}

#' Hierarchical list for K=9 rescue test.
hier_k9 <- function(K = 9L) {
  list(
    coarse_mask      = c(rep(1L, 3L), rep(0L, K - 3L)),
    min_cell_n       = 30L,
    mode             = 0L,
    outer_tol        = 1e-4,
    outer_iterations = 100L
  )
}
