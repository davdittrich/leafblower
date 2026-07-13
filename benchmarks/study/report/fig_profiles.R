# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# fig_profiles.R -- family-stratified Dolan-More performance profiles over
# the S3C parametric instance family (if_* problems), read-only from
# benchmarks/study/results/metrics.parquet (cost = wall_time_s_median; no
# metric recompute per S5/I3).
#
# Method: Dolan & More 2002 (Math. Prog. 91(2):201-213,
# doi:10.1007/s101070100263) -- performance ratio r_{p,s} = t_{p,s} /
# min_s' t_{p,s'}; a DNF/failed run gets the sentinel r = infinity, realised
# here as r_M = 1 + max(finite ratio) *within the stratum*. The profile
# rho_s(tau) = |{p : r_{p,s} <= tau}| / |P| is exactly the empirical CDF of
# solver s's ratio sample -> plotted via an ECDF evaluated on a log2(tau)
# grid (ggplot2::stat_ecdf's own per-panel grid is not reused here because
# the bootstrap ribbon below needs the SAME evaluation grid as the point
# estimate; see the tau_grid comment).
#
# Gould & Scott 2016 (ACM TOMS 43(2):15, doi:10.1145/2950048) S7: profiles
# are problem-set-sensitive -> compute *within* each objective-family
# stratum (kl / chi2 / logit / minimax, the `family` column / `strata` in
# _data.R), never pooled, and attach a bootstrap-over-PROBLEMS (not reps)
# pointwise 2.5/97.5% CI ribbon (B >= 1000) so a reader can see how much of
# a stratum's profile is supported by how few problems.
#
# `family` (metrics.parquet) is a per-SOLVER attribute -- the divergence
# that solver's objective targets -- constant across every row of that
# solver (confirmed by inspection 2026-07-13: `n_families == 1` for all 39
# solvers appearing on if_* problems). A stratum's solver set is therefore
# fixed by `family`, not chosen per problem.
#
# |P| per stratum = number of distinct if_* problems with >=1 logged
# attempt (by any solver) under that family -- NOT a fixed 30 for every
# panel: e.g. the `logit` stratum only ever ran 1 of the 30 if_* problems
# in this snapshot, a validity-floor breach (see VALIDITY_FLOOR below) that
# this script surfaces loudly (console + caption) rather than hiding it by
# reusing the full 30-problem denominator.
#
# DNF definition: status %in% DNF_STATUS, OR missing weight/timing columns
# ("no real weights", defensive -- does not trigger on the current
# snapshot), OR a (solver, problem) cell entirely absent from the family's
# full solver x problem grid (a solver that never attempted a problem is
# still a Dolan-More failure on it, not silently excluded from the
# denominator).
#
# Build-variant discipline (DESIGN.md S5): `build_series(variant =
# "headline")` restricts leafblower to build == "portable" (confirmed:
# leafblower's Python bindings are the portable builds and its R bindings
# are the -march=native builds in this snapshot) so each leafblower method
# enters the headline profile exactly ONCE, never mixed with its
# build == "native" rows -- mixing would double-count leafblower and
# contaminate the per-problem min-cost normalizer and rho_s(1). A second,
# clearly-labelled native-delta profile (leafblower-only, build ==
# "native", competitors dropped) is written separately, faceted by family
# in one combined PDF (facet_wrap only -- patchwork is not installed).

.candidate_report_dirs <- c("benchmarks/study/report", "study/report", ".")
.this_dir <- Find(function(d) file.exists(file.path(d, "_data.R")), .candidate_report_dirs)
if (is.null(.this_dir)) stop("fig_profiles.R: could not locate report/_data.R from working directory: ", getwd())

source(file.path(.this_dir, "_data.R"))

FIG_DIR <- file.path(.this_dir, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# --- constants -------------------------------------------------------------

DNF_STATUS <- c("dnf", "error", "infeasible", "bound_violation", "budget", "stall", "no_conv")
PROFILE_B <- 1000L       # bootstrap-over-problems replicates (Gould & Scott S7: B >= 1000)
PROFILE_SEED <- 20260713L
VALIDITY_FLOOR <- 20L    # |P| below this is flagged, not silently rendered as if reliable
TAU_GRID_N <- 300L       # log2(tau) evaluation grid for both the point estimate and the CI ribbon

# --- load + restrict to the S3C instance family (if_*) ---------------------

metrics <- read_metrics()
metrics_if <- metrics[grepl("^if_", metrics$problem), , drop = FALSE]
stopifnot("no if_* (S3C instance family) rows found in metrics.parquet" = nrow(metrics_if) > 0)

flag_dnf <- function(df) {
  df$dnf <- df$status %in% DNF_STATUS |
    is.na(df$wall_time_s_median) | is.na(df$wmin) | is.na(df$wmax) | is.na(df$weight_kl)
  df
}

headline <- flag_dnf(build_series(metrics_if, variant = "headline"))
native_delta <- flag_dnf(build_series(metrics_if, variant = "native_delta"))

n_p_all <- length(unique(metrics_if$problem))
message(sprintf("[fig_profiles] S3C instance family: |P|_total (all if_* problems) = %d", n_p_all))

# --- core: build one family's ratio grid + ECDF profile + bootstrap CI -----
#
# Returns a long-format data.frame (family, solver, tier, stale, log2tau,
# rho, ci_lo, ci_hi) plus metadata (n_p, n_solver, r_M, solvers).

build_family_profile <- function(df, family_id) {
  fam <- df[df$family == family_id, , drop = FALSE]
  stopifnot(
    "fig_profiles: >1 row per (solver,problem) in family stratum (thread>1? build not collapsed?)" =
      !any(duplicated(fam[c("solver", "problem")]))
  )
  solvers <- sort(unique(fam$solver))
  problems <- sort(unique(fam$problem))
  n_solver <- length(solvers)
  n_p <- length(problems)
  stopifnot(
    "family has no solvers" = n_solver > 0,
    "family has no problems" = n_p > 0
  )

  # Full solver x problem grid: a missing cell means the solver never
  # attempted that problem -- still a DNF, not an omission from |P|.
  grid <- expand.grid(solver = solvers, problem = problems, stringsAsFactors = FALSE)
  grid <- dplyr::left_join(
    grid,
    fam[, c("solver", "problem", "wall_time_s_median", "dnf")],
    by = c("solver", "problem")
  )
  grid$dnf[is.na(grid$dnf)] <- TRUE
  grid$cost <- ifelse(grid$dnf, NA_real_, grid$wall_time_s_median)

  non_dnf <- grid[!grid$dnf, , drop = FALSE]
  if (nrow(non_dnf) > 0) {
    min_cost <- stats::aggregate(cost ~ problem, data = non_dnf, FUN = min)
    names(min_cost)[2] <- "min_cost"
  } else {
    min_cost <- data.frame(problem = character(0), min_cost = numeric(0))
  }
  grid <- dplyr::left_join(grid, min_cost, by = "problem")

  # NA (dnf, or every solver dnf on that problem, i.e. min_cost undefined)
  # is replaced by the stratum-wide sentinel r_M below -- no 0/0 or Inf
  # ever reaches the plotted ratio.
  grid$ratio <- grid$cost / grid$min_cost
  finite_ratio <- grid$ratio[is.finite(grid$ratio)]
  r_M <- if (length(finite_ratio) > 0) 1 + max(finite_ratio) else 2
  grid$ratio_sentinel <- ifelse(is.finite(grid$ratio), grid$ratio, r_M)

  mat <- matrix(NA_real_, nrow = n_p, ncol = n_solver, dimnames = list(problems, solvers))
  mat[cbind(match(grid$problem, problems), match(grid$solver, solvers))] <- grid$ratio_sentinel
  stopifnot("ratio_sentinel grid has unfilled cells" = !anyNA(mat))

  # Common log2(tau) evaluation grid, shared by the point estimate and the
  # bootstrap ribbon (a fixed, evenly-spaced grid -- not the raw data
  # break points -- keeps bootstrap memory O(B * TAU_GRID_N * n_solver)
  # instead of O(B * n_p * n_solver * n_solver)). rho(tau) at each grid
  # point is an EXACT count via findInterval (no interpolation).
  log2tau_grid <- seq(0, log2(r_M), length.out = TAU_GRID_N)
  tau_grid <- 2^log2tau_grid

  ecdf_at <- function(v) findInterval(tau_grid, sort(v)) / n_p

  rho_point <- apply(mat, 2, ecdf_at)  # TAU_GRID_N x n_solver

  set.seed(PROFILE_SEED)
  boot_arr <- array(NA_real_, dim = c(PROFILE_B, TAU_GRID_N, n_solver))
  for (b in seq_len(PROFILE_B)) {
    idx <- sample.int(n_p, n_p, replace = TRUE)
    boot_arr[b, , ] <- apply(mat[idx, , drop = FALSE], 2, ecdf_at)
  }
  ci_lo <- apply(boot_arr, c(2, 3), stats::quantile, probs = 0.025, names = FALSE, type = 7)
  ci_hi <- apply(boot_arr, c(2, 3), stats::quantile, probs = 0.975, names = FALSE, type = 7)

  curve_df <- do.call(rbind, lapply(seq_len(n_solver), function(j) {
    data.frame(
      family = family_id,
      solver = solvers[j],
      log2tau = log2tau_grid,
      rho = rho_point[, j],
      ci_lo = ci_lo[, j],
      ci_hi = ci_hi[, j],
      stringsAsFactors = FALSE
    )
  }))
  curve_df$tier <- tier_of(curve_df$solver)
  curve_df$stale <- stale_of(curve_df$solver)

  attempted <- vapply(solvers, function(s) sum(!grid$dnf[grid$solver == s]), integer(1))
  below_floor <- n_p < VALIDITY_FLOOR
  message(sprintf(
    "[fig_profiles] family=%-8s |P|=%2d solvers=%2d B=%d r_M=%.3f%s",
    family_id, n_p, n_solver, PROFILE_B, r_M,
    if (below_floor) "  *** VALIDITY FLOOR BREACH (|P| < 20, Gould-Scott S7) ***" else ""
  ))

  list(
    curve_df = curve_df,
    n_p = n_p,
    n_solver = n_solver,
    r_M = r_M,
    below_floor = below_floor,
    attempted = attempted
  )
}

# --- shared caption text -----------------------------------------------------

.caption_common <- function(prof, variant_label) {
  floor_note <- if (prof$below_floor) {
    sprintf(
      "\n*** VALIDITY-FLOOR WARNING: |P| = %d < %d (Gould & Scott 2016 S7) -- this stratum's profile/CI is NOT reliable, most solvers were only attempted on a small fraction of the S3C instance family. ***",
      prof$n_p, VALIDITY_FLOOR
    )
  } else ""
  paste0(
    "Dolan & More (2002) performance profile: r = wall_time_s_median / min-over-solvers, DNF/absent-cell -> r = r_M sentinel (1 + max finite ratio in stratum). ",
    "x = log2(tau); rho_s(tau) = ECDF of r for solver s. ",
    "Family-stratified per Gould & Scott (2016) S7 (profiles are problem-set-sensitive) -- computed within this objective-family stratum only, NOT pooled across strata. ",
    sprintf("|P| = %d if_* problems, %d solvers, B = %d bootstrap-over-PROBLEMS replicates (pointwise 2.5/97.5%% ribbon). ", prof$n_p, prof$n_solver, PROFILE_B),
    variant_label, " ",
    "Colour = competitor-discipline tier; dashed line = [!] solver labelled 'historical baseline' (stale maintenance, DESIGN.md S7).",
    floor_note
  )
}

# --- render: one PDF per family (headline, portable leafblower) ------------

headline_profiles <- list()
for (fam in strata) {
  prof <- build_family_profile(headline, fam)
  headline_profiles[[fam]] <- prof
  cd <- prof$curve_df

  p <- ggplot2::ggplot(cd, ggplot2::aes(x = log2tau, y = rho, group = solver, colour = tier)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lo, ymax = ci_hi, fill = tier), alpha = 0.12, colour = NA) +
    ggplot2::geom_step(ggplot2::aes(linetype = stale), linewidth = 0.5) +
    ggplot2::scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    theme_lbw() +
    ggplot2::labs(
      title = sprintf("Dolan-Moré performance profile -- %s family (S3C instance suite, headline/portable)", fam),
      subtitle = sprintf("|P| = %d if_* problems, %d solvers, B = %d bootstrap replicates", prof$n_p, prof$n_solver, PROFILE_B),
      x = expression(log[2](tau)),
      y = expression(rho[s](tau)),
      colour = "tier",
      fill = "tier",
      caption = .caption_common(prof, "Headline build variant: leafblower restricted to build == \"portable\" via build_series() so each method enters once.")
    )

  out_path <- file.path(FIG_DIR, paste0("profile_", fam, ".pdf"))
  ggplot2::ggsave(out_path, p, width = 10, height = 7)
  stopifnot(
    "profile figure was not written" = file.exists(out_path),
    "profile figure is empty" = file.info(out_path)$size > 0
  )
  message("wrote ", out_path)
}

# --- render: optional native-delta profile (leafblower-only, build == "native") --
# One combined PDF, faceted by family (facet_wrap only -- patchwork not
# installed): each panel only has the small number of leafblower solvers
# that shipped a native build in this snapshot (competitors dropped
# entirely by build_series(variant = "native_delta") -- there is no native
# competitor variant to compare against).

native_curves <- lapply(strata, function(fam) {
  prof <- build_family_profile(native_delta, fam)
  prof$curve_df
})
native_df <- do.call(rbind, native_curves)

p_native <- ggplot2::ggplot(native_df, ggplot2::aes(x = log2tau, y = rho, group = solver, colour = tier)) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ci_lo, ymax = ci_hi, fill = tier), alpha = 0.12, colour = NA) +
  ggplot2::geom_step(ggplot2::aes(linetype = stale), linewidth = 0.5) +
  ggplot2::scale_linetype_manual(values = c(`FALSE` = "solid", `TRUE` = "dashed"), guide = "none") +
  ggplot2::coord_cartesian(ylim = c(0, 1)) +
  ggplot2::facet_wrap(~family, scales = "free_x") +
  theme_lbw() +
  ggplot2::labs(
    title = "Dolan-Moré performance profile -- native-build delta (leafblower-only, S3C instance suite)",
    subtitle = sprintf("B = %d bootstrap-over-PROBLEMS replicates per family panel; competitors excluded (no native competitor variant)", PROFILE_B),
    x = expression(log[2](tau)),
    y = expression(rho[s](tau)),
    colour = "tier",
    fill = "tier",
    caption = paste0(
      "Dolan & More (2002) performance profile computed leafblower-only, within build == \"native\" ",
      "(build_series(variant = \"native_delta\")) -- NOT a headline comparison against competitors. ",
      "Family-stratified per Gould & Scott (2016) S7, facet_wrap only (patchwork not installed). ",
      "DNF/absent-cell -> r_M sentinel per family panel; [!] dashed = 'historical baseline' stale solver."
    )
  )

out_path_native <- file.path(FIG_DIR, "profile_native_delta.pdf")
ggplot2::ggsave(out_path_native, p_native, width = 11, height = 8)
stopifnot(
  "native-delta profile figure was not written" = file.exists(out_path_native),
  "native-delta profile figure is empty" = file.info(out_path_native)$size > 0
)
message("wrote ", out_path_native)

# --- optional derived ratio parquet (report/_derived/, gitignored) --------

source(file.path(.this_dir, "_derived.R"))

derived_headline <- do.call(rbind, lapply(headline_profiles, function(p) p$curve_df))
write_derived(derived_headline, "profile_ratios_headline")
write_derived(native_df, "profile_ratios_native_delta")
