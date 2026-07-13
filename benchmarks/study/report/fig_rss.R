# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# fig_rss.R — peak-RSS (RQ2/RQ3) figure + table, read-only from
# benchmarks/study/results/metrics.parquet$peak_rss_bytes.
#
# peak_rss_bytes is ALREADY baseline-subtracted at data-load time (§5) and is
# a max-over-reps high-water mark per (solver, problem, build) row — it is
# NOT a repeated-measures sample, so no median/mean±sd is computed here, only
# the (already-final) per-row value. `build_series(variant = "headline")` is
# used to pick exactly one row per (solver, problem) (the CRAN/PyPI-default
# leafblower build vs. each competitor's single build), per the DESIGN.md §5
# chokepoint documented in _data.R — this avoids double-counting leafblower's
# portable+native rows against single-build competitors.

.candidate_report_dirs <- c("benchmarks/study/report", "study/report", ".")
.this_dir <- Find(function(d) file.exists(file.path(d, "_data.R")), .candidate_report_dirs)
if (is.null(.this_dir)) stop("fig_rss.R: could not locate report/_data.R from working directory: ", getwd())

source(file.path(.this_dir, "_data.R"))

FIG_DIR <- file.path(.this_dir, "figures")
TBL_DIR <- file.path(.this_dir, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# --- language lookup (registry.json `arm` field: "R" or "python") ---------
# Read-only lookup of the ground-truth per-solver language, not a guess from
# the solver id string (competitor ids like "svy"/"anesrake" don't encode
# language in the id).
REGISTRY_PATH <- file.path(dirname(RESULTS_DIR), "registry.json")
if (!file.exists(REGISTRY_PATH)) stop("fig_rss.R: registry.json not found: ", REGISTRY_PATH)
.registry <- jsonlite::fromJSON(REGISTRY_PATH, simplifyVector = FALSE)
.solver_lang <- data.frame(
  solver   = names(.registry$solvers),
  language = vapply(.registry$solvers, function(x) if (!is.null(x$arm)) x$arm else NA_character_, character(1)),
  stringsAsFactors = FALSE
)
.solver_lang$language <- ifelse(tolower(.solver_lang$language) == "r", "R",
                          ifelse(tolower(.solver_lang$language) == "python", "Python", .solver_lang$language))

metrics <- read_metrics()
rss <- build_series(metrics, variant = "headline") |>
  dplyr::left_join(.solver_lang, by = "solver") |>
  dplyr::mutate(
    package        = solver,
    peak_rss_mib   = peak_rss_bytes / 1024^2,
    tier           = tier_of(solver)
  )

stopifnot(
  "peak_rss_bytes missing after join" = !anyNA(rss$peak_rss_bytes),
  "language lookup left unmatched solvers" = !anyNA(rss$language)
)

# --- figure: per-(solver, problem) peak RSS, log10 y (bytes span ~120 MiB
# to ~36 GiB across problem sizes n=1e3..1.58e6) ---------------------------

p <- ggplot2::ggplot(rss, ggplot2::aes(x = stats::reorder(solver, peak_rss_mib), y = peak_rss_mib, colour = language)) +
  ggplot2::geom_point(size = 1.6) +
  ggplot2::scale_y_log10(name = "peak RSS (MiB, baseline-subtracted, log10)") +
  ggplot2::facet_wrap(~problem, scales = "free_x") +
  ggplot2::coord_flip() +
  theme_lbw() +
  ggplot2::labs(
    title = "Peak RSS by solver and problem",
    subtitle = "Max-over-reps high-water mark; baseline already subtracted at data-load time",
    x = "solver",
    colour = "language",
    caption = paste(
      "peak_rss_bytes is a per-(solver,problem,build) max over repetitions,",
      "not a repeated-measures sample -- no median/mean+-sd is reported.",
      "Baseline RSS was already subtracted upstream (metrics.parquet); this",
      "figure does not re-subtract it. One row per (solver,problem) via",
      "build_series(variant=\"headline\").",
      sep = "\n"
    )
  )

fig_path <- file.path(FIG_DIR, "rss.pdf")
ggplot2::ggsave(fig_path, p, width = 11, height = 8)
stopifnot(
  "rss figure was not written" = file.exists(fig_path),
  "rss figure is empty" = file.info(fig_path)$size > 0
)
message("wrote ", fig_path)

# --- table: per-(solver, problem) peak RSS, plus a per-(language, package)
# rollup (max over problems, since these are already high-water marks and
# taking a further max preserves the "high-water mark" semantics without
# introducing a median/mean) --------------------------------------------

tbl_by_solver_problem <- rss |>
  dplyr::arrange(problem, dplyr::desc(peak_rss_mib)) |>
  dplyr::select(problem, solver, language, tier, peak_rss_mib, peak_rss_bytes)

tbl_by_language_package <- rss |>
  dplyr::group_by(language, package) |>
  dplyr::summarise(
    peak_rss_mib_max = max(peak_rss_mib),
    n_problems        = dplyr::n(),
    .groups = "drop"
  ) |>
  dplyr::arrange(dplyr::desc(peak_rss_mib_max))

utils::write.csv(tbl_by_solver_problem, file.path(TBL_DIR, "rss_by_solver_problem.csv"), row.names = FALSE)
utils::write.csv(tbl_by_language_package, file.path(TBL_DIR, "rss_by_language_package.csv"), row.names = FALSE)

for (f in c("rss_by_solver_problem.csv", "rss_by_language_package.csv")) {
  fp <- file.path(TBL_DIR, f)
  stopifnot(
    "rss table was not written" = file.exists(fp),
    "rss table is empty" = file.info(fp)$size > 0
  )
  message("wrote ", fp)
}
