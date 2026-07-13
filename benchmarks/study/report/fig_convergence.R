# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# fig_convergence.R — per-iteration convergence curves for the 8 instrumented
# solvers, read directly from benchmarks/study/results/trajectories/*.csv
# (read-only; no metric recompute per §5/I3).
#
# Trajectory CSV schema (confirmed by inspection, 2026-07-13):
#   filename  <solver>__<problem>__<lang>.csv   (lang in {py, r})
#   columns   iter,<natural-metric>,marginal_kl   (sparse checkpoints:
#             1,2,3,5,10,20,50,100,200,500 — the exact checkpoint set varies
#             slightly per run; do not assume a fixed grid, just plot what's
#             present). 2-column (iter,<metric>) files are still accepted
#             (marginal_kl absent -> NA) for forward/backward tolerance.
#
#   Column 2 (the "natural" metric) is each solver's OWN internal
#   stopping-criterion residual/gap proxy, NOT uniform across solvers --
#   confirmed by reading every CSV header on disk:
#     errRp        — chebyshev, greenkhorn, oris, oris_soft, raking, sinkhorn
#     dual_gap     — newton_kl
#     chi2         — greg
#   These are mixed units and NOT comparable across panels with a different
#   metric name -- see the headline figures' caption.
#
#   Column 3, `marginal_kl` (present in all 64 on-disk files, including
#   newton_kl, all finite), is the metric harvest() reports for the KL
#   family: a common, cross-solver-comparable convergence axis. It is
#   monotone but can show floating-point cancellation noise (~-1e-14) right
#   at machine-precision convergence; those points are clamped to a small
#   positive display floor for the log10 axis only (see below) -- no value
#   is recomputed.
#
# Covers exactly the 8 instrumented solvers with trajectory files on disk;
# no other solver is plotted (they have no per-iteration data).

.candidate_report_dirs <- c("benchmarks/study/report", "study/report", ".")
.this_dir <- Find(function(d) file.exists(file.path(d, "_data.R")), .candidate_report_dirs)
if (is.null(.this_dir)) stop("fig_convergence.R: could not locate report/_data.R from working directory: ", getwd())

source(file.path(.this_dir, "_data.R"))

TRAJ_DIR <- file.path(RESULTS_DIR, "trajectories")
if (!dir.exists(TRAJ_DIR)) stop("fig_convergence.R: trajectories dir not found: ", TRAJ_DIR)

FIG_DIR <- file.path(.this_dir, "figures")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# Display-only floor for the marginal_kl log10 axis: FP-cancellation noise in
# the reconstructed KL sits around -1e-14 at machine-precision convergence;
# clamping to machine epsilon keeps those points visible near the axis floor
# without fabricating a value or recomputing the metric.
MARGINAL_KL_DISPLAY_FLOOR <- .Machine$double.eps

# --- load + parse every trajectory CSV -----------------------------------

traj_files <- list.files(TRAJ_DIR, pattern = "\\.csv$", full.names = TRUE)
if (length(traj_files) == 0) stop("fig_convergence.R: no trajectory CSVs found in ", TRAJ_DIR)

# Filename grammar: <solver>__<problem>__<lang>.csv (solver/problem may
# themselves contain single underscores, so split only on the DOUBLE
# underscore separator: last segment = lang, first segment = solver, the
# (possibly multi-segment) middle = problem).
.parse_traj_name <- function(path) {
  base <- sub("\\.csv$", "", basename(path))
  parts <- strsplit(base, "__", fixed = TRUE)[[1]]
  if (length(parts) < 3) stop("fig_convergence.R: unexpected trajectory filename: ", basename(path))
  lang <- parts[length(parts)]
  solver <- parts[1]
  problem <- paste(parts[2:(length(parts) - 1)], collapse = "__")
  list(solver = solver, problem = problem, lang = lang)
}

# Each CSV has iter, <natural-metric> and (since 2026-07-13) a 3rd
# `marginal_kl` column; 2-column files (no marginal_kl) are still tolerated.
# Normalise to `metric_value`/`metric_name` (natural) plus `marginal_kl`
# (NA if the column is absent).
.read_one_traj <- function(f) {
  meta <- .parse_traj_name(f)
  d <- utils::read.csv(f)
  stopifnot(
    "trajectory CSV must have 2 (iter, <metric>) or 3 (iter, <metric>, marginal_kl) columns" =
      ncol(d) %in% c(2, 3),
    "trajectory CSV's first column must be 'iter'" = names(d)[1] == "iter"
  )
  metric_name <- names(d)[2]
  if (ncol(d) == 3) {
    stopifnot(
      "trajectory CSV's 3rd column must be 'marginal_kl'" = names(d)[3] == "marginal_kl"
    )
    marginal_kl <- as.double(d[[3]])
  } else {
    marginal_kl <- NA_real_
  }
  data.frame(
    iter = as.integer(d[[1]]),
    metric_value = as.double(d[[2]]),
    metric_name = metric_name,
    marginal_kl = marginal_kl,
    solver = meta$solver,
    problem = meta$problem,
    lang = meta$lang,
    stringsAsFactors = FALSE
  )
}

traj <- do.call(rbind, lapply(traj_files, .read_one_traj))

# The 8 instrumented solvers (all trajectory files on disk belong to these;
# assert rather than silently drop unexpected ids).
INSTRUMENTED_SOLVERS <- c(
  "chebyshev", "greenkhorn", "greg", "newton_kl",
  "oris", "oris_soft", "raking", "sinkhorn"
)
stopifnot(
  "unexpected solver id in trajectories/ (not in the 8 instrumented solvers)" =
    all(traj$solver %in% INSTRUMENTED_SOLVERS)
)
# One natural metric name per solver, uniformly across all its files
# (checked, not assumed): fail loudly if that invariant ever breaks.
.metric_per_solver <- unique(traj[, c("solver", "metric_name")])
stopifnot(
  "a solver has more than one distinct natural trajectory metric column across its files" =
    !any(duplicated(.metric_per_solver$solver))
)

traj$lang <- factor(traj$lang, levels = c("r", "py"), labels = c("R", "Python"))
traj$facet_label <- paste0(traj$solver, " [", traj$metric_name, "]")

problems <- sort(unique(traj$problem))

# --- (a) HEADLINE figures: natural default metric per solver, one faceted --
# figure per problem. facet_wrap only (patchwork not installed) -- one panel
# per solver (each panel's own metric named in its strip), colour = solver,
# linetype = lang. Honest labelling: mixed units, NOT comparable across
# panels with a different metric name.

for (prob in problems) {
  d <- traj[traj$problem == prob, ]
  d <- d[order(d$solver, d$lang, d$iter), ]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = iter, y = metric_value, colour = solver, linetype = lang)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~facet_label, scales = "free_y") +
    theme_lbw() +
    ggplot2::labs(
      title = paste0("Convergence trajectories — ", prob),
      subtitle = "8 instrumented solvers; sparse checkpoints; log10 y-axis",
      x = "iteration",
      y = "solver-internal residual/gap proxy (log10) -- see facet strip for metric name",
      colour = "solver",
      linetype = "language",
      caption = paste(
        "Facet strips name each solver's own internal stopping-criterion",
        "quantity (errRp / dual_gap / chi2) -- these are NOT the marginal-KL",
        "divergence used for the cross-solver comparable figure below, and",
        "are not comparable across panels with different metric names.",
        "Checkpoints are sparse (not every iteration was logged).",
        sep = "\n"
      )
    )

  out_path <- file.path(FIG_DIR, paste0("convergence_", prob, ".pdf"))
  ggplot2::ggsave(out_path, p, width = 10, height = 7)
  stopifnot(
    "convergence figure was not written" = file.exists(out_path),
    "convergence figure is empty" = file.info(out_path)$size > 0
  )
  message("wrote ", out_path)
}

# --- (b) marginal_kl COMPARABLE figures: common cross-solver axis, one -----
# faceted figure per problem (facet = solver only, since the metric is the
# same in every panel). Includes newton_kl (its marginal_kl is finite, unlike
# its natural dual_gap which lives on a different scale).

stopifnot(
  "marginal_kl column missing for one or more trajectory rows -- cannot build the comparable figure" =
    all(!is.na(traj$marginal_kl))
)

traj$marginal_kl_display <- pmax(traj$marginal_kl, MARGINAL_KL_DISPLAY_FLOOR)

for (prob in problems) {
  d <- traj[traj$problem == prob, ]
  d <- d[order(d$solver, d$lang, d$iter), ]

  p <- ggplot2::ggplot(d, ggplot2::aes(x = iter, y = marginal_kl_display, colour = solver, linetype = lang)) +
    ggplot2::geom_line() +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_y_log10() +
    ggplot2::facet_wrap(~solver, scales = "free_y") +
    theme_lbw() +
    ggplot2::labs(
      title = paste0("Marginal-KL convergence — ", prob),
      subtitle = "8 instrumented solvers; sparse checkpoints; log10 y-axis",
      x = "iteration",
      y = "marginal_kl (log10)",
      colour = "solver",
      linetype = "language",
      caption = paste(
        "marginal_kl is the common convergence axis (the metric harvest()",
        "reports for the KL family): monotone and comparable across ALL",
        "solvers, including newton_kl. Values are clamped to a",
        sprintf("%.2e", MARGINAL_KL_DISPLAY_FLOOR),
        "display floor for the log10 axis only, to show FP-cancellation",
        "noise (~-1e-14) at machine-precision convergence without",
        "recomputing the metric. Checkpoints are sparse.",
        sep = "\n"
      )
    )

  out_path <- file.path(FIG_DIR, paste0("convergence_marginal_kl_", prob, ".pdf"))
  ggplot2::ggsave(out_path, p, width = 10, height = 7)
  stopifnot(
    "marginal_kl convergence figure was not written" = file.exists(out_path),
    "marginal_kl convergence figure is empty" = file.info(out_path)$size > 0
  )
  message("wrote ", out_path)
}
