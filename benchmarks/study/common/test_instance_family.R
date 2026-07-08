#!/usr/bin/env Rscript
# test_instance_family.R -- tests for benchmarks/study/common/instance_family.R
# (WU-3, leafblower-2ouc.4).
#
# Usage: Rscript benchmarks/study/common/test_instance_family.R [roundtrip_out_path]
# When roundtrip_out_path is given, writes a compact JSON dump of a few small
# instances' generated category matrices + targets for the companion Python
# test to independently regenerate and diff against (R<->Py exact-parity
# check, DoD item 4).
suppressPackageStartupMessages({
  library(jsonlite)
})

here <- function(...) file.path("benchmarks", "study", ...)
source(here("common", "problem_io.R"))
source(here("common", "instance_family.R"))

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

cat("== build_instance_grid() ==\n")
grid <- build_instance_grid()
check("30 <= n_instances <= 100", length(grid) >= 30L && length(grid) <= 100L)
cat(sprintf("  (n_instances = %d)\n", length(grid)))

n_vals <- vapply(grid, function(x) x$n, numeric(1))
K_vals <- vapply(grid, function(x) x$K, numeric(1))
card_vals <- vapply(grid, function(x) x$card, character(1))
skew_vals <- vapply(grid, function(x) x$skew, character(1))
cond_vals <- vapply(grid, function(x) x$cond, character(1))
infeas_vals <- vapply(grid, function(x) x$infeas, character(1))

check("n axis covers all 4 frozen levels", setequal(unique(n_vals), .IF_N_LEVELS))
check("K axis covers all 3 frozen levels", setequal(unique(K_vals), .IF_K_LEVELS))
check("card axis covers low/medium/high", setequal(unique(card_vals), names(.IF_CARD_LEVELS)))
check("skew axis covers none/moderate/extreme", setequal(unique(skew_vals), names(.IF_SKEW_LEVELS)))
check("cond axis covers well/moderate/ill", setequal(unique(cond_vals), names(.IF_COND_LEVELS)))
check("infeas axis covers loose/moderate/tight", setequal(unique(infeas_vals), names(.IF_INFEAS_MAX)))
check("only a handful (<=4) of instances touch n=1.58M", sum(n_vals == 1580000L) <= 4L)
check("no duplicate instance ids", !any(duplicated(vapply(grid, .if_instance_id, character(1)))))

cat("== determinism: .if_generate_categories is repeatable ==\n")
g1 <- .if_generate_categories(200L, 4L, 8L, 1, 0.5, 0L)
g2 <- .if_generate_categories(200L, 4L, 8L, 1, 0.5, 0L)
check("identical matrices across repeated calls", identical(g1$mat, g2$mat))
g3 <- .if_generate_categories(200L, 4L, 8L, 1, 0.5, 1L)
check("different seed => different matrix", !identical(g1$mat, g3$mat))
check("category values within [1,C]", all(g1$mat >= 1L) && all(g1$mat <= 8L))

cat("== generate_instance_family_specs() writes schema-conformant specs ==\n")
tmp_dir <- file.path(tempdir(), "instance_family_spec_test")
res <- generate_instance_family_specs(tmp_dir)
check("n written ids == grid length", length(res$ids) == length(grid))
sample_ids <- res$ids[c(1L, round(length(res$ids) / 2), length(res$ids))]
for (id in sample_ids) {
  spec <- jsonlite::fromJSON(file.path(res$out_dir, paste0(id, ".json")), simplifyVector = FALSE)
  required <- c("id", "data_ref", "margins", "bounds", "tol", "objective_families", "K")
  check(paste0(id, ": has all required fields"), all(required %in% names(spec)))
  check(paste0(id, ": data_ref starts with gen:instance_family"),
        grepl("^gen:instance_family\\?", spec$data_ref))
  margins <- unlist(spec$margins)
  check(paste0(id, ": length(margins) == K"), length(margins) == spec$K)
  check(paste0(id, ": targets keys == margins"), setequal(names(spec$targets), margins))
  for (m in margins) {
    check(paste0(id, ": targets$", m, " sums to 1"),
          near(sum(unlist(spec$targets[[m]])), 1))
  }
  check(paste0(id, ": bounds$min == 0"), near(spec$bounds$min, 0))
}

cat("== gen: resolution through the WU-2 loader (load_problem_spec) ==\n")
install_gen_resolver()
small_ids <- res$ids[n_vals[seq_along(res$ids)] <= 10000L]
for (id in small_ids[seq_len(min(3L, length(small_ids)))]) {
  path <- file.path(res$out_dir, paste0(id, ".json"))
  problem <- load_problem_spec(path)
  spec_raw <- jsonlite::fromJSON(path, simplifyVector = FALSE)
  check(paste0(id, ": loaded n matches spec n query param"),
        nrow(problem$data) == as.integer(sub(".*[?&]n=([0-9]+).*", "\\1", spec_raw$data_ref)))
  check(paste0(id, ": K margin columns present"), all(problem$margins %in% names(problem$data)))
  check(paste0(id, ": margin columns are factors"), all(vapply(problem$margins, function(m) is.factor(problem$data[[m]]), logical(1))))
  check(paste0(id, ": design_weights default to all-ones"), all(problem$design_weights == 1))
}

cat("== unrecognized gen: recipe still falls through to WU-2's not-implemented guard ==\n")
gen_err <- tryCatch({
  .pio_resolve_data_ref("gen:toy_recipe", NULL)
  NULL
}, error = function(e) conditionMessage(e))
check("gen:toy_recipe (unregistered) still raises WU-3-scope error after install_gen_resolver()",
      is.character(gen_err) && grepl("WU-3", gen_err))

if (length(commandArgs(trailingOnly = TRUE)) >= 1L) {
  out_path <- commandArgs(trailingOnly = TRUE)[1]
  roundtrip_specs <- list(
    list(n = 60L, K = 3L, C = 5L, s = 1, rho = 0.5, seed = 0L),
    list(n = 40L, K = 2L, C = 4L, s = 3, rho = 0.9, seed = 7L)
  )
  dumps <- lapply(roundtrip_specs, function(p) {
    g <- .if_generate_categories(p$n, p$K, p$C, p$s, p$rho, p$seed)
    list(
      n = p$n, K = p$K, C = p$C, s = p$s, rho = p$rho, seed = p$seed,
      mat = matrix(as.integer(g$mat), nrow = p$n, ncol = p$K),
      probs = as.numeric(g$zipf$probs)
    )
  })
  # digits = 17 -- max decimal significant digits needed to round-trip an
  # IEEE-754 double exactly (DBL_DIG+... margin); digits=12 (prior value)
  # truncated Zipf probs enough to break the EXACT parity assertion below.
  writeLines(jsonlite::toJSON(dumps, auto_unbox = TRUE, digits = 17), out_path)
  cat(sprintf("\nWrote round-trip dump to %s\n", out_path))
}

cat(sprintf("\nRESULT: %d failure(s)\n", failures))
if (failures > 0L) quit(status = 1L)
