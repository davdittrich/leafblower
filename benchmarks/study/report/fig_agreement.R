# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# fig_agreement.R — RQ5 weight-vector agreement heatmaps, STRICTLY-CONVEX
# families only (kl, chi2, logit — unique optimum). Minimax/L-inf is
# EXCLUDED here by construction (DESIGN.md §6 "Blocker G": the L-inf LP
# optimum lies on a face, not a unique vertex, so two correct solvers return
# different weight vectors at identical achieved L-inf — Pearson would
# falsely report disagreement); minimax's objective-value agreement TABLE
# (rel_gap) lives in tbl_quality.R::minimax_objval, not here.
#
# One tile-heatmap PDF per family: report/figures/agreement_<family>.pdf.
# Axes = solver x solver (all pairwise combinations present in
# agreement.parquet — analysis/aggregate_metrics.py:259 computes every
# unordered pair per problem, not solely vs. a single reference row). Fill =
# median Pearson correlation of raw weight vectors across problem instances
# (median, not mean — matches the report-wide median+quantile convention).
# `degenerate` cells (zero-variance/NaN pearson on >=1 side, per
# common/metrics.py:agreement()) are masked to a DISTINCT grey — never
# scored as 0 or 1 agreement.
#
# Read-only on results/agreement.parquet + results/metrics.parquet (the
# latter only for its already-recomputed W/n columns, to build the RQ5
# "mean(w)" normaliser below) — no metric recompute.

.candidate_report_dirs <- c("benchmarks/study/report", "study/report", ".")
.this_dir <- Find(function(d) file.exists(file.path(d, "_fig_common.R")), .candidate_report_dirs)
if (is.null(.this_dir)) stop("fig_agreement.R: could not locate report/_fig_common.R from working directory: ", getwd())
source(file.path(.this_dir, "_fig_common.R"))
suppressPackageStartupMessages(library(dplyr))

STRICTLY_CONVEX_FAMILIES <- c("kl", "chi2", "logit")

# RQ5 pre-registered "tight" thresholds (DESIGN.md §1/§6).
RQ5_PEARSON_MIN  <- 0.999
RQ5_SPEARMAN_MIN <- 0.999
RQ5_RELDIFF_MAX  <- 1e-3

agreement <- read_agreement()
metrics   <- read_metrics()

# Headline build discipline (per _data.R build_series() chokepoint comment):
# an agreement row carries build_a/build_b for TWO independent solvers, so
# build_series()'s single-`build`-column contract doesn't apply verbatim —
# this is its two-sided equivalent: keep only rows where BOTH sides are a
# competitor's single build ("na") or leafblower's CRAN/PyPI-default
# ("portable"), so a leafblower solver's "native" (-O3 -march=native) cell
# never mixes into the headline agreement pool.
.is_headline_build <- function(build_col) build_col %in% c("na", "portable")

wv <- agreement |>
  filter(mode == "weight_vector", family %in% STRICTLY_CONVEX_FAMILIES) |>
  filter(.is_headline_build(build_a), .is_headline_build(build_b))

stopifnot(
  "no strictly-convex weight_vector rows found after headline-build filter" = nrow(wv) > 0
)

# mean(w) per (solver, problem, build) = W/n — already-recomputed columns
# (no recompute), NOT assumed == 1: the project's Sigma(w)=n normalisation
# is enforced for leafblower but several competitor packages (ebal, sbw,
# nonprobsvy, ...) do not exit at that normalisation (verified against
# metrics.parquet 2026-07-13: mean_w ranges 1e-8..182 across solvers) — so
# this join is required to make the RQ5 "max|Delta w|/mean(w)" relative
# threshold meaningful rather than silently assuming mean(w)==1.
mean_w <- metrics |> transmute(solver, problem, build, mean_w = W / n)

wv <- wv |>
  left_join(mean_w |> rename(solver_a = solver, build_a = build, mean_w_a = mean_w),
            by = c("solver_a", "problem", "build_a")) |>
  left_join(mean_w |> rename(solver_b = solver, build_b = build, mean_w_b = mean_w),
            by = c("solver_b", "problem", "build_b")) |>
  mutate(rel_max_abs_diff = max_abs_diff / ((mean_w_a + mean_w_b) / 2))

stopifnot(
  "mean_w join left unmatched (solver,problem,build) rows -- metrics.parquet coverage gap" =
    !anyNA(wv$mean_w_a) && !anyNA(wv$mean_w_b)
)

# Per-(family, solver-pair) aggregate across problem instances. A row's
# `degenerate` flag is `!(isfinite(pearson) & isfinite(spearman))`
# (common/metrics.py:agreement()) — broader than "pearson is NaN": verified
# against agreement.parquet (2026-07-13), 413 of 691 kl-family degenerate
# rows carry a spuriously finite near-zero `pearson` (a 0/0 floating-point
# artifact of a zero-variance side; e.g. pearson=4.47e-16) while `spearman`
# is the one that resolves to NaN. Including those spurious values in
# `median_na(pearson)` would silently pollute a partially-degenerate pair's
# median with noise, so `!degenerate` filters ALL THREE aggregated columns
# BEFORE taking the median (same subset of well-defined problem-instances
# per column) — pearson_med is NA_real_ **only** when frac_degenerate == 1
# (every compared problem was degenerate for that pair; median_na(numeric(0))
# == NA), which is exactly the "degenerate cell" this script masks.
cell_agg <- wv |>
  group_by(family, solver_a, solver_b) |>
  summarise(
    n_problems      = dplyr::n(),
    frac_degenerate = mean(degenerate),
    pearson_med     = median_na(pearson[!degenerate]),
    spearman_med    = median_na(spearman[!degenerate]),
    rel_diff_med    = median_na(rel_max_abs_diff[!degenerate]),
    .groups = "drop"
  ) |>
  mutate(
    meets_rq5 = !is.na(pearson_med) & !is.na(spearman_med) & !is.na(rel_diff_med) &
      pearson_med >= RQ5_PEARSON_MIN & spearman_med >= RQ5_SPEARMAN_MIN &
      rel_diff_med <= RQ5_RELDIFF_MAX
  )

# Symmetrise for a full square tile grid: agreement.parquet stores each
# unordered pair once (itertools.combinations(..., 2)); mirror (b,a) so it
# renders identically to (a,b) — this is a pure display-layer duplication,
# not a second independent observation (no double-counting in cell_agg's
# own median, which is computed before this mirror step).
cell_sym <- dplyr::bind_rows(
  cell_agg,
  cell_agg |> dplyr::rename(solver_a = solver_b, solver_b = solver_a)
)

#' Build one family's Pearson-agreement tile heatmap.
plot_family_heatmap <- function(family_name) {
  d <- cell_sym |> dplyr::filter(family == family_name)
  stopifnot("plot_family_heatmap: no rows for this family" = nrow(d) > 0)

  solvs <- sort(unique(c(d$solver_a, d$solver_b)))
  lbl <- stats::setNames(stale_label(solvs), solvs)
  d <- d |> mutate(
    solver_a = factor(solver_a, levels = solvs),
    solver_b = factor(solver_b, levels = rev(solvs))
  )
  n_pairs_headline <- nrow(cell_agg |> dplyr::filter(family == family_name))
  max_problems <- max(d$n_problems)

  ggplot2::ggplot(d, ggplot2::aes(solver_a, solver_b, fill = pearson_med)) +
    ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(
      ggplot2::aes(label = ifelse(meets_rq5, "✓", "")),
      colour = "black", size = 3
    ) +
    ggplot2::scale_fill_gradient2(
      name = "median\nPearson r",
      low = "firebrick", mid = "white", high = "steelblue",
      midpoint = 0, limits = c(-1, 1),
      na.value = "grey60"
    ) +
    ggplot2::scale_x_discrete(labels = lbl) +
    ggplot2::scale_y_discrete(labels = lbl) +
    theme_lbw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = paste0("Weight-vector agreement — ", family_name, " family (strictly-convex)"),
      subtitle = paste0(
        length(solvs), " solvers, ", n_pairs_headline, " pairs, up to ", max_problems,
        " problem instance(s) per pair (headline build only)"
      ),
      caption = paste(
        "Strictly-convex families only (kl/chi2/logit -- unique optimum); minimax is EXCLUDED",
        "(judged on objective-value agreement instead -- see tbl_quality.R::minimax_objval, DESIGN.md Blocker G).",
        "Grey tile = degenerate/undefined (zero-variance weight vector on >=1 side, ALL compared",
        "problem instances) -- NOT scored as 0 or 1 agreement.",
        "✓ = RQ5 tight threshold met: median Pearson>=0.999 AND Spearman>=0.999 AND",
        "max|Delta w|/mean(w)<=1e-3 (max_abs_diff normalised by the pair's mean weight scale).",
        "⚠ = stale/unmaintained package (historical baseline snapshot, not actively updated).",
        "Headline build selection only (leafblower: portable; competitors: single build) --",
        "see _data.R::build_series().",
        sep = "\n"
      )
    )
}

for (family_name in STRICTLY_CONVEX_FAMILIES) {
  p <- plot_family_heatmap(family_name)
  n_solv <- length(unique(c(
    cell_sym$solver_a[cell_sym$family == family_name],
    cell_sym$solver_b[cell_sym$family == family_name]
  )))
  side <- max(6, 0.4 * n_solv + 3)
  fig_path <- file.path(FIG_DIR, paste0("agreement_", family_name, ".pdf"))
  ggplot2::ggsave(fig_path, p, width = side + 2, height = side, device = grDevices::cairo_pdf, limitsize = FALSE)
  stopifnot(
    "agreement figure was not written" = file.exists(fig_path),
    "agreement figure is empty" = file.info(fig_path)$size > 0
  )
  message("wrote ", fig_path)
}
