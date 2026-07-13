# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# fig_work_precision.R — RQ1 work-precision figures: one log-log scatter
# PDF per objective family (kl, chi2, logit, minimax) of achieved
# native-divergence error (x) vs. median end-to-end wall time (y), one
# point per (solver, problem), colored by tier_of(solver), shaped by
# build-variant (leafblower headline = portable build, native = separate
# "tuned-build delta" overlay per build_series()). The empirical Pareto
# envelope (sort-on-error, running-min-time — see _fig_common.R::
# pareto_front()) is drawn as a black step line, never eyeballed.
#
# Family -> error column (inspected against results/metrics.parquet on
# 2026-07-13; native_div_kind's actual distinct values are chi2_dist,
# logit_dist, margin_linf, unknown, weight_kl — NOT the "chi2"/"logit"
# guessed in the WU brief):
#   kl      -> marg_kl_mean                 (native_div differs by orders
#                                             of magnitude for kl; NOT used)
#   chi2    -> native_div  (native_div_kind == "chi2_dist")
#   logit   -> native_div  (native_div_kind == "logit_dist")
#   minimax -> margin_linf (NEUTRAL axis only — Chebyshev/minimax solvers
#                            are home-field advantaged on this metric;
#                            see caption caveat)
# The 8 family == "dispatch" rows (native_div_kind == "unknown") are
# outside `strata` (kl/chi2/logit/minimax) and are excluded, matching
# _data.R's `strata` definition.
#
# Read-only on results/metrics.parquet via _data.R::read_metrics(); no
# metric recompute (only column selection + a display-only log floor, see
# _fig_common.R::floor_for_log()). Renders headless: `Rscript
# fig_work_precision.R` from repo root or from report/.

.candidate_report_dirs <- c("benchmarks/study/report", "study/report", ".")
.this_dir <- Find(function(d) file.exists(file.path(d, "_fig_common.R")), .candidate_report_dirs)
if (is.null(.this_dir)) stop("fig_work_precision.R: could not locate report/_fig_common.R from working directory: ", getwd())
source(file.path(.this_dir, "_fig_common.R"))

#' Compute the per-family native-divergence error column + a human-readable
#' provenance label for axis titles/captions.
error_for_family <- function(df) {
  stopifnot(
    "chi2 family rows must have native_div_kind == 'chi2_dist' (re-inspect metrics.parquet before editing this mapping)" =
      all(df$native_div_kind[df$family == "chi2"] == "chi2_dist"),
    "logit family rows must have native_div_kind %in% c('logit_dist','logit')" =
      all(df$native_div_kind[df$family == "logit"] %in% c("logit_dist", "logit"))
  )
  df |>
    mutate(
      err = case_when(
        family == "kl"      ~ marg_kl_mean,
        family == "chi2"    ~ native_div,
        family == "logit"   ~ native_div,
        family == "minimax" ~ margin_linf,
        TRUE ~ NA_real_
      ),
      err_label = case_when(
        family == "kl"      ~ "Marginal KL divergence (marg_kl_mean)",
        family == "chi2"    ~ "Native chi-squared divergence (native_div, chi2_dist)",
        family == "logit"   ~ "Native logit divergence (native_div, logit_dist)",
        family == "minimax" ~ "L-infinity margin (margin_linf) -- NEUTRAL axis",
        TRUE ~ NA_character_
      )
    )
}

metrics_all <- read_metrics() |>
  build_series(variant = "all") |>
  filter(family %in% strata) |>
  error_for_family() |>
  mutate(
    err = ave(err, family, FUN = floor_for_log),
    tier = tier_of(solver),
    stale = stale_of(solver)
  )

for (fam in strata) {

  fam_df <- metrics_all |> filter(family == fam, build_role %in% c("competitor", "headline", "delta"))
  if (nrow(fam_df) == 0) {
    warning("fig_work_precision.R: no rows for family '", fam, "' — skipping")
    next
  }

  envelope <- pareto_front(fam_df)
  err_label <- unique(fam_df$err_label)[1]

  p <- ggplot(fam_df, aes(x = err, y = wall_time_s_median)) +
    geom_point(aes(color = tier, shape = build_role), size = 2.4, alpha = 0.85) +
    geom_step(
      data = envelope, aes(x = err, y = envelope_time),
      inherit.aes = FALSE, direction = "hv", color = "black", linewidth = 0.6
    ) +
    scale_x_log10(labels = scales::label_log()) +
    scale_y_log10(labels = scales::label_log()) +
    scale_color_manual(values = TIER_COLORS, name = "Tier") +
    scale_shape_manual(values = BUILD_SHAPES, name = "Build variant") +
    labs(
      title = paste0("Work precision -- ", fam, " family"),
      subtitle = "Black step line = empirical Pareto envelope (sort-on-error, running-min-time)",
      x = err_label,
      y = "Median end-to-end wall time (s, log scale)",
      caption = paste0(
        "Accuracy-normalised; leafblower build variant: portable = headline, ",
        "native = tuned-build delta (open triangle). ",
        if (any(fam_df$stale)) {
          paste0("⚠ historical baseline (stale packages, labelled with † elsewhere): ",
                 paste(sort(unique(fam_df$solver[fam_df$stale])), collapse = ", "), ". ")
        } else "",
        if (fam == "minimax") {
          "minimax family: margin_linf used as a NEUTRAL error axis only -- Chebyshev/minimax solvers are home-field advantaged on this metric; do not read it as a cross-family accuracy ranking."
        } else ""
      )
    ) +
    theme_lbw()

  out_path <- file.path(FIG_DIR, paste0("work_precision_", fam, ".pdf"))
  ggsave(out_path, p, width = 7.5, height = 5.5, device = cairo_pdf)
  stopifnot(
    "work-precision PDF must render non-empty" = file.exists(out_path) && file.info(out_path)$size > 0
  )
  message("wrote ", out_path, " (", nrow(fam_df), " points, ", nrow(envelope), " on Pareto envelope)")
}
