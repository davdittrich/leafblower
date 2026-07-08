#!/usr/bin/env Rscript
# test_problem_io.R — loader smoke test for benchmarks/study/common/problem_io.R
#
# Exercises all four data_ref origins reachable from R:
#   inline  -> spec/toy_inline.json
#   pkg:    -> spec/canonical_survey_apistrat.json (survey::apistrat)
#   file:   -> spec/stepstone_unbounded.json, spec/stepstone_bounded.json
#   gen:    -> asserted to raise the documented WU-3-not-implemented error
#
# Usage: Rscript benchmarks/study/common/test_problem_io.R [roundtrip_out_path]
# When roundtrip_out_path is given, writes a compact JSON summary of the two
# cross-language-portable specs (toy_inline, stepstone_unbounded) for the
# companion Python test to diff against (R<->Py round-trip check, WU-2 DoD).
suppressPackageStartupMessages({
  library(jsonlite)
})

here <- function(...) file.path("benchmarks", "study", ...)
source(here("common", "problem_io.R"))

failures <- 0L
check <- function(desc, cond) {
  if (isTRUE(cond)) {
    cat(sprintf("  PASS: %s\n", desc))
  } else {
    cat(sprintf("  FAIL: %s\n", desc))
    failures <<- failures + 1L
  }
}
near <- function(a, b, tol = 1e-9) abs(a - b) <= tol

cat("== toy_inline (data_ref='inline') ==\n")
toy <- load_problem_spec(here("spec", "toy_inline.json"))
check("id", identical(toy$id, "toy_inline"))
check("n rows", nrow(toy$data) == 4L)
check("K", toy$K == 1L)
check("margins", identical(toy$margins, "grp"))
check("targets sum to 1", near(sum(toy$targets$grp), 1))
check("targets values", near(toy$targets$grp[["A"]], 0.5) && near(toy$targets$grp[["B"]], 0.5))
check("design_weights (inline array)", identical(toy$design_weights, c(1, 1, 2, 2)))
check("bounds min", near(toy$bounds$min, 0))
check("bounds max", near(toy$bounds$max, 10))
check("margin column is factor", is.factor(toy$data$grp))

cat("== canonical_survey_apistrat (data_ref='pkg:survey::apistrat') ==\n")
api <- load_problem_spec(here("spec", "canonical_survey_apistrat.json"))
check("id", identical(api$id, "canonical_survey_apistrat"))
check("n rows", nrow(api$data) == 200L)
check("K", api$K == 1L)
check("margins", identical(api$margins, "stype"))
pop_total <- 4421 + 755 + 1018
check("targets sum to 1", near(sum(api$targets$stype), 1))
check("targets E proportion", near(api$targets$stype[["E"]], 4421 / pop_total))
check("targets H proportion", near(api$targets$stype[["H"]], 755 / pop_total))
check("targets M proportion", near(api$targets$stype[["M"]], 1018 / pop_total))
check("design_weights length", length(api$design_weights) == 200L)
check("design_weights matches pw sum", near(sum(api$design_weights), 6193.99995804, tol = 1e-4))
check("bounds max unbounded", is.infinite(api$bounds$max))

cat("== stepstone_unbounded (data_ref='file:...parquet') ==\n")
step_u <- load_problem_spec(here("spec", "stepstone_unbounded.json"))
check("id", identical(step_u$id, "stepstone_unbounded"))
check("n rows", nrow(step_u$data) == 1582732L)
check("K", step_u$K == 9L)
check("margins count", length(step_u$margins) == 9L)
check("design_weights all ones", all(step_u$design_weights == 1))
check("design_weights length", length(step_u$design_weights) == nrow(step_u$data))
check("bounds max unbounded", is.infinite(step_u$bounds$max))
for (nm in step_u$margins) {
  check(paste0("targets '", nm, "' sum to 1"), near(sum(step_u$targets[[nm]]), 1))
}

cat("== stepstone_bounded (data_ref='file:...parquet', bounds.max=5) ==\n")
step_b <- load_problem_spec(here("spec", "stepstone_bounded.json"))
check("id", identical(step_b$id, "stepstone_bounded"))
check("n rows", nrow(step_b$data) == 1582732L)
check("bounds max = 5", near(step_b$bounds$max, 5))

cat("== gen: origin (WU-3 not-yet-implemented guard) ==\n")
gen_err <- tryCatch({
  .pio_resolve_data_ref("gen:toy_recipe", NULL)
  NULL
}, error = function(e) conditionMessage(e))
check("gen: raises WU-3-scope error", is.character(gen_err) && grepl("WU-3", gen_err))

if (!is.null(step_u$data) && length(commandArgs(trailingOnly = TRUE)) >= 1L) {
  out_path <- commandArgs(trailingOnly = TRUE)[1]
  # jsonlite pitfalls guarded here: (1) auto_unbox=TRUE unboxes length-1
  # atomic vectors to scalars, so a K=1 `margins` character(1) would
  # serialize as a bare string, not a 1-element array -- wrap with I() to
  # force array output regardless of length. (2) toJSON drops names on a
  # plain named numeric vector (renders {"A":.5,"B":.5} as [.5,.5]) -- convert
  # each target vector via as.list() first so names survive as JSON object
  # keys. (3) toJSON renders Inf as the string "Inf"; map to NA (-> JSON
  # null) so the schema's own null-means-unbounded convention round-trips.
  .fmt_bound <- function(x) if (is.infinite(x)) NA_real_ else x
  summary_obj <- list(
    toy_inline = list(
      id = toy$id, n = nrow(toy$data), K = toy$K, margins = I(toy$margins),
      targets = lapply(toy$targets, as.list),
      design_weights_sum = sum(toy$design_weights),
      bounds = list(min = .fmt_bound(toy$bounds$min), max = .fmt_bound(toy$bounds$max))
    ),
    stepstone_unbounded = list(
      id = step_u$id, n = nrow(step_u$data), K = step_u$K,
      margins = I(sort(step_u$margins)),
      targets = lapply(step_u$targets[sort(names(step_u$targets))], as.list),
      design_weights_sum = sum(step_u$design_weights),
      bounds = list(min = .fmt_bound(step_u$bounds$min), max = .fmt_bound(step_u$bounds$max))
    )
  )
  writeLines(jsonlite::toJSON(summary_obj, auto_unbox = TRUE, digits = 12), out_path)
  cat(sprintf("\nWrote round-trip summary to %s\n", out_path))
}

cat(sprintf("\n%s: %d failure(s)\n", if (failures == 0L) "RESULT" else "RESULT", failures))
if (failures > 0L) quit(status = 1L)
