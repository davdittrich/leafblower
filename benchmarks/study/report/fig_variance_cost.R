# STUDY-BRANCH-ONLY-DO-NOT-MERGE
#
# fig_variance_cost.R — RQ6 downstream variance/replicate-weight cost
# figure + table, read-only from
# benchmarks/study/results/variance_cost.parquet.
#
# variance_cost.parquet is a per-(solver, problem) summary of a B=200
# Rao-Wu rescaled-bootstrap COLD re-solve loop (each replicate re-solves the
# calibration dual from scratch, mirroring survey::calibrate.svyrep.design's
# reset-per-replicate contract) — NOT a repeated-measures sample of a single
# solve, so this script reports the already-final per-(solver,problem)
# summary columns (median/min/max of per-replicate wall time, the resulting
# total B-replicate cost, convergence counts, and the inflation_factor
# cross-check) rather than re-deriving them.

.candidate_report_dirs <- c("benchmarks/study/report", "study/report", ".")
.this_dir <- Find(function(d) file.exists(file.path(d, "_data.R")), .candidate_report_dirs)
if (is.null(.this_dir)) stop("fig_variance_cost.R: could not locate report/_data.R from working directory: ", getwd())

source(file.path(.this_dir, "_data.R"))
source(file.path(.this_dir, "_fig_common.R"))

FIG_DIR <- file.path(.this_dir, "figures")
TBL_DIR <- file.path(.this_dir, "tables")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

#' Read the RQ6 variance-estimation-cost summary (B=200 cold replicate
#' re-solves per (solver, problem)). Read-only — no metric recompute, same
#' contract as read_metrics()/read_runs() in _data.R.
read_variance_cost <- function(path = file.path(RESULTS_DIR, "variance_cost.parquet")) {
  arrow::read_parquet(path)
}

vc <- read_variance_cost() |>
  dplyr::mutate(tier = tier_of(solver))

stopifnot(
  "fig_variance_cost: >1 row per (solver,problem) -- grain must be one row per (solver,problem)" =
    !any(duplicated(vc[c("solver", "problem")]))
)

# --- figure: total_variance_cost (B=200 cold re-solves) per solver, faceted
# by problem, log10 y (leafblower rows visually distinct via tier colour,
# same TIER_COLORS/order_solver_factor convention as the other fig_*.R) ----

vc_ord <- order_solver_factor(vc, tier_col = "tier", time_col = "total_variance_cost")

p <- ggplot2::ggplot(vc_ord, ggplot2::aes(x = solver, y = total_variance_cost, fill = tier)) +
  ggplot2::geom_col() +
  ggplot2::scale_y_log10(name = "total variance-estimation cost (s, B=200 cold re-solves, log10)") +
  ggplot2::scale_fill_manual(values = TIER_COLORS) +
  ggplot2::facet_wrap(~problem, scales = "free_y") +
  ggplot2::coord_flip() +
  theme_lbw() +
  ggplot2::labs(
    title = "Downstream variance-estimation cost by solver and problem",
    subtitle = "Rao-Wu rescaled-bootstrap, B=200 COLD replicate re-solves (dual reset per replicate)",
    x = "solver",
    fill = "tier",
    caption = paste(
      "total_variance_cost = B * per_replicate_wall_median; one row per",
      "(solver, problem). Cold re-solve mirrors",
      "survey::calibrate.svyrep.design's per-replicate dual reset -- a",
      "solver cheap to fit once may rank differently once re-solved B",
      "times for variance estimation.",
      sep = "\n"
    )
  )

fig_path <- file.path(FIG_DIR, "variance_cost.pdf")
ggplot2::ggsave(fig_path, p, width = 11, height = 8)
stopifnot(
  "variance_cost figure was not written" = file.exists(fig_path),
  "variance_cost figure is empty" = file.info(fig_path)$size > 0
)
message("wrote ", fig_path)

# --- table: per (solver, problem) cost + convergence + cross-check --------

tbl <- vc |>
  dplyr::arrange(problem, dplyr::desc(total_variance_cost)) |>
  dplyr::mutate(n_converged_of_B = paste0(n_converged, "/", B)) |>
  dplyr::select(
    problem, solver, tier,
    total_variance_cost, per_replicate_wall_median,
    n_converged_of_B, inflation_factor
  )

tbl_path <- file.path(TBL_DIR, "variance_cost.csv")
utils::write.csv(tbl, tbl_path, row.names = FALSE)
stopifnot(
  "variance_cost table was not written" = file.exists(tbl_path),
  "variance_cost table is empty" = file.info(tbl_path)$size > 0
)
message("wrote ", tbl_path)
