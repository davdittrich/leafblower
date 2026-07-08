#!/usr/bin/env Rscript
# benchmarks/study/common/produce_ref_rows.R
#
# Production driver (ticket leafblower-2ouc.5, WU-4): solves+stores the
# independent convex reference anchor for every currently-resolvable
# canonical problem spec, for each family that spec declares in its
# `objective_families` (intersected with REF_FAMILIES). NOT a test -- this
# writes the actual pseudo-solver rows consumed by RQ5.
#
# Currently-resolvable canonical specs (see benchmarks/study/spec/*.json):
#   toy_inline.json                -- resolvable in R and Python (data_ref=inline)
#   canonical_survey_apistrat.json -- resolvable in R ONLY (data_ref=pkg:survey::
#                                      apistrat, an R-only package; problem_io.py
#                                      has no Python resolver for pkg: by design)
# stepstone_bounded/unbounded.json (n ~ 1.58M) are DELIBERATELY EXCLUDED --
# REF_MAX_N in ref_convex.R refuses them with an explicit error rather than
# faking a ~1e-12 anchor at a scale where that is infeasible (DESIGN.md
# Section 6). instance-family (`gen:` origin) specs are WU-3 scope, not yet
# implemented -- none exist to produce rows for yet.
#
# Run: Rscript benchmarks/study/common/produce_ref_rows.R

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
.script_dir <- if (length(.script)) dirname(.script) else "benchmarks/study/common"
source(file.path(.script_dir, "ref_convex.R"))
source(file.path(.script_dir, "problem_io.R"))

repo_root <- normalizePath(file.path(.script_dir, "..", "..", ".."))
spec_dir <- file.path(repo_root, "benchmarks", "study", "spec")
out_dir <- file.path(repo_root, "benchmarks", "study", "results")

canonical_specs <- c("toy_inline.json", "canonical_survey_apistrat.json")

for (spec_file in canonical_specs) {
  spec_path <- file.path(spec_dir, spec_file)
  problem <- load_problem_spec(spec_path)
  families <- intersect(problem$objective_families, REF_FAMILIES)
  cat(sprintf("== %s (n=%d): families=%s ==\n", problem$id, nrow(problem$data),
              paste(families, collapse = ",")))
  for (family in families) {
    res <- tryCatch(solve_ref(problem, family), error = function(e) {
      cat(sprintf("  [SKIP] family=%s: %s\n", family, conditionMessage(e)))
      NULL
    })
    if (is.null(res)) next
    path <- store_ref(problem, family, res, out_dir)
    cat(sprintf("  [OK] family=%-8s mode=%-16s status=%-10s -> %s\n",
                family, res$mode, res$solver_status, path))
  }
}

cat("\n== stepstone (n ~ 1.58M) -- documented as having NO anchor, not faked ==\n")
for (spec_file in c("stepstone_unbounded.json", "stepstone_bounded.json")) {
  spec_path <- file.path(spec_dir, spec_file)
  ss_problem <- load_problem_spec(spec_path)  # genuine load -- real n, not a stand-in
  err <- tryCatch({ solve_ref(ss_problem, "kl"); NULL }, error = function(e) conditionMessage(e))
  cat(sprintf("  [%s] n=%d guard message: %s\n", ss_problem$id, nrow(ss_problem$data), err))
}
