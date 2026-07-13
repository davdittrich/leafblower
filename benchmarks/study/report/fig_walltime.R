# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# fig_walltime.R — RQ1 end-to-end wall-time figure: one PDF per objective
# family (kl, chi2, logit, minimax), a single panel per family with one row
# per solver (median time as a point, [median-of-mins, median-of-maxes] as
# a 5/95-proxy whisker per wall_time_s_min/wall_time_s_max — NO mean/sd
# anywhere), tier-colored, leafblower headline (portable build) vs. native
# ("tuned-build delta") distinguished by shape per build_series().
#
# Per-solver aggregation is across that family's problems, using
# _data.R::median_na() (pure median, matching this study's "median only"
# invariant) — NOT a re-implementation of build filtering (build_series()
# is the sole chokepoint, called once below).
#
# Read-only on results/metrics.parquet via _data.R::read_metrics(); no
# metric recompute (grouping/summarising of already-recomputed columns
# only). Renders headless: `Rscript fig_walltime.R` from repo root or from
# report/.

.candidate_report_dirs <- c("benchmarks/study/report", "study/report", ".")
.this_dir <- Find(function(d) file.exists(file.path(d, "_fig_common.R")), .candidate_report_dirs)
if (is.null(.this_dir)) stop("fig_walltime.R: could not locate report/_fig_common.R from working directory: ", getwd())
source(file.path(.this_dir, "_fig_common.R"))

metrics_all <- read_metrics() |>
  build_series(variant = "all") |>
  filter(family %in% strata, build_role %in% c("competitor", "headline", "delta")) |>
  mutate(tier = tier_of(solver))

for (fam in strata) {

  fam_raw <- metrics_all |> filter(family == fam)
  if (nrow(fam_raw) == 0) {
    warning("fig_walltime.R: no rows for family '", fam, "' — skipping")
    next
  }

  stopifnot(
    "fig_walltime: >1 row per (solver,problem,build_role) in fam_raw (thread>1?) -- across-problems median would fold threads" =
      !any(duplicated(fam_raw[c("solver", "problem", "build_role")]))
  )

  fam_agg <- fam_raw |>
    group_by(solver, tier, build_role) |>
    summarise(
      wall_time_s_median = median_na(wall_time_s_median),
      wall_time_p5        = median_na(wall_time_s_min),
      wall_time_p95       = median_na(wall_time_s_max),
      .groups = "drop"
    ) |>
    as.data.frame()

  fam_agg <- order_solver_factor(fam_agg)
  fam_agg$solver_label <- factor(stale_label(as.character(fam_agg$solver)),
                                  levels = stale_label(levels(fam_agg$solver)))

  n_rows <- nrow(fam_agg)
  stale_note <- if (any(stale_of(fam_agg$solver))) {
    paste0("⚠ historical baseline (stale packages): ",
           paste(sort(unique(fam_agg$solver[stale_of(fam_agg$solver)])), collapse = ", "), ". ")
  } else ""
  minimax_note <- if (fam == "minimax") {
    "minimax family: recall margin_linf (used elsewhere as this family's error axis) is a NEUTRAL metric only -- Chebyshev/minimax solvers are home-field advantaged on it. "
  } else ""

  p <- ggplot(fam_agg, aes(y = solver_label)) +
    geom_segment(aes(x = wall_time_p5, xend = wall_time_p95, yend = solver_label), linewidth = 3, alpha = 0.25, color = "grey40") +
    geom_point(aes(x = wall_time_s_median, color = tier, shape = build_role), size = 2.6) +
    scale_x_log10(labels = scales::label_log()) +
    scale_color_manual(values = TIER_COLORS, name = "Tier") +
    scale_shape_manual(values = BUILD_SHAPES, name = "Build variant") +
    labs(
      title = paste0("End-to-end wall time -- ", fam, " family"),
      subtitle = "Point = median across problems; band = median-of-min to median-of-max (5/95 proxy), NOT mean±sd",
      x = "Wall time (s, log scale)",
      y = NULL,
      caption = paste0(
        "Accuracy-normalised study design (all runs at a common nominal tolerance); ",
        "leafblower build variant: portable = headline, native = tuned-build delta (open triangle). ",
        stale_note, minimax_note
      )
    ) +
    theme_lbw()

  out_path <- file.path(FIG_DIR, paste0("walltime_", fam, ".pdf"))
  ggsave(out_path, p, width = 7.5, height = max(3.5, 1.2 + 0.28 * n_rows), device = cairo_pdf, limitsize = FALSE)
  stopifnot(
    "wall-time PDF must render non-empty" = file.exists(out_path) && file.info(out_path)$size > 0
  )
  message("wrote ", out_path, " (", n_rows, " solvers)")
}
