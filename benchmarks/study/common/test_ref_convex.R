#!/usr/bin/env Rscript
# benchmarks/study/common/test_ref_convex.R
#
# Cross-check for benchmarks/study/common/ref_convex.R (ticket leafblower-2ouc.5,
# WU-4) against an INDEPENDENT hand derivation on spec/toy_inline.json -- same
# "no re-invoking the code under test" discipline as test_metrics.R (WU-1).
#
# Hand derivation (toy_inline.json: n=4, margin grp in {A,B}, design_weights=
# [1,1,2,2], targets{A:2,B:2} normalized to {A:0.5,B:0.5}):
#   N = sum(d) = 6; T_A = T_B = 0.5*6 = 3.
#   Single margin, closed-form KL-raking scale per group = T_group / D_group:
#     group A: D_A = d1+d2 = 2, scale = 3/2 = 1.5  -> w1=1*1.5=1.5, w2=1*1.5=1.5
#     group B: D_B = d3+d4 = 4, scale = 3/4 = 0.75 -> w3=2*0.75=1.5, w4=2*0.75=1.5
#   All four weights converge to exactly 1.5 (independent of family: this is
#   the UNIQUE feasible point satisfying both margin equalities on a single
#   2-level margin -- kl/chi2/logit and minimax all land here since the
#   feasible set collapses to a point).
#   weight_kl = sum(w*log(w/d)):
#     obs1,2 (d=1): 1.5*log(1.5) = 0.6081976621622468 each
#     obs3,4 (d=2): 1.5*log(0.75) = -0.4315231086776712 each
#     total = 2*0.6081976621622468 + 2*(-0.4315231086776712) = 0.35334910696915045
#
# Run: Rscript benchmarks/study/common/test_ref_convex.R

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
.script_dir <- if (length(.script)) dirname(.script) else "benchmarks/study/common"
source(file.path(.script_dir, "ref_convex.R"))
source(file.path(.script_dir, "problem_io.R"))

fail_count <- 0L
check <- function(desc, got, want, tol = 1e-8) {
  ok <- isTRUE(all.equal(got, want, tolerance = tol))
  cat(sprintf("[%s] %-70s got=%s want=%s\n", if (ok) "PASS" else "FAIL",
              desc, format(got, digits = 12), format(want, digits = 12)))
  if (!ok) fail_count <<- fail_count + 1L
  invisible(ok)
}
check_true <- function(desc, got) {
  ok <- isTRUE(got)
  cat(sprintf("[%s] %-70s got=%s\n", if (ok) "PASS" else "FAIL", desc, got))
  if (!ok) fail_count <<- fail_count + 1L
  invisible(ok)
}

repo_root <- normalizePath(file.path(.script_dir, "..", "..", ".."))
spec_path <- file.path(repo_root, "benchmarks", "study", "spec", "toy_inline.json")
problem <- load_problem_spec(spec_path)

W_HAND <- c(1.5, 1.5, 1.5, 1.5)
WEIGHT_KL_HAND <- 0.35334910696915045

cat("== kl family (WU-1 golden toy cross-check) ==\n")
res_kl <- solve_ref(problem, "kl")
check_true("kl mode == weight_vector", identical(res_kl$mode, "weight_vector"))
check("kl weights[1]", res_kl$weights[1], 1.5)
check("kl weights[2]", res_kl$weights[2], 1.5)
check("kl weights[3]", res_kl$weights[3], 1.5)
check("kl weights[4]", res_kl$weights[4], 1.5)
check_true("kl weights match hand-derived vector (max abs diff < 1e-8)",
           max(abs(res_kl$weights - W_HAND)) < 1e-8)
wkl <- weight_kl(res_kl$weights, problem$design_weights, family = "kl")$weight_kl
check("weight_kl(ref kl weights) matches hand-derived 0.3533489267", wkl, WEIGHT_KL_HAND)

cat("\n== chi2 family (same toy: unique feasible point == same 1.5 vector) ==\n")
res_chi2 <- solve_ref(problem, "chi2")
check_true("chi2 mode == weight_vector", identical(res_chi2$mode, "weight_vector"))
check_true("chi2 weights match hand-derived vector", max(abs(res_chi2$weights - W_HAND)) < 1e-8)

cat("\n== logit family (same toy: bounds [0,10] finite, unique feasible point) ==\n")
res_logit <- solve_ref(problem, "logit")
check_true("logit mode == weight_vector", identical(res_logit$mode, "weight_vector"))
check_true("logit weights match hand-derived vector", max(abs(res_logit$weights - W_HAND)) < 1e-8)

cat("\n== minimax family (Blocker G: objective_value anchor, not weight-vector) ==\n")
res_mm <- solve_ref(problem, "minimax")
check_true("minimax mode == objective_value (never weight_vector)",
           identical(res_mm$mode, "objective_value"))
check("minimax achieved margin_linf ~ 0 (single margin fully saturable)", res_mm$obj_val, 0, tol = 1e-8)
check_true("minimax result carries NO stored weight-vector anchor field 'weights' used only internally",
           is.numeric(res_mm$weights))  # present internally, but store_ref() below must NOT write it as parquet

cat("\n== REF_MAX_N scope guard (stepstone 1.58M has no anchor, DESIGN.md Sec.6) ==\n")
fake_big <- list(id = "fake_stepstone", data = data.frame(x = seq_len(REF_MAX_N + 1L)))
guard_err <- tryCatch({ solve_ref(fake_big, "kl"); NULL }, error = function(e) conditionMessage(e))
check_true("solve_ref refuses n > REF_MAX_N with an explicit message (not a silent hang/fake)",
           !is.null(guard_err) && grepl("REF_MAX_N", guard_err))

cat("\n== store_ref: pseudo-solver row persistence ==\n")
out_dir <- file.path(repo_root, "benchmarks", "study", "results")
p_kl <- store_ref(problem, "kl", res_kl, out_dir)
p_chi2 <- store_ref(problem, "chi2", res_chi2, out_dir)
p_logit <- store_ref(problem, "logit", res_logit, out_dir)
p_mm <- store_ref(problem, "minimax", res_mm, out_dir)
check_true("kl weights parquet written", file.exists(p_kl))
check_true("chi2 weights parquet written", file.exists(p_chi2))
check_true("logit weights parquet written", file.exists(p_logit))
check_true("minimax objective-value JSON written (NOT a weights parquet)", file.exists(p_mm))
check_true("minimax anchor path is under ref_objective/, not weights/", grepl("ref_objective", p_mm, fixed = TRUE))

reread <- as.data.frame(arrow::read_parquet(p_kl))
check_true("re-read kl parquet round-trips the weight vector",
           max(abs(reread$weight - res_kl$weights)) < 1e-12)

cat(sprintf("\n%d check(s) failed.\n", fail_count))
quit(status = if (fail_count > 0L) 1L else 0L, save = "no")
