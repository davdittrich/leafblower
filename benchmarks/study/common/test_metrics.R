#!/usr/bin/env Rscript
# benchmarks/study/common/test_metrics.R
#
# Hand-computed golden for benchmarks/study/common/metrics.R (ticket leafblower-2ouc.2).
# Golden values are computed by an INDEPENDENT scratch derivation (fractions
# for exact cases; log() as an arithmetic oracle for KL terms, applied to
# hand-derived exact p/T vectors) -- not by re-invoking the code under test.
#
# Also validates the THREE legacy metric impls (verbatim function bodies,
# copied for reference -- the source files are NOT modified) against the same
# toy: the base (non-divergent) case must agree with the new impl on
# marg_kl/max_err/L1 (undisputed formulas); the starved-category and d_i!=1
# cases are expected to REPRODUCE the known bugs (silently dropped divergent
# term; hardcoded d_i=1) -- this is the golden actively demonstrating why the
# refactor is required, per DESIGN.md Section 6.
#
# Run: Rscript benchmarks/study/common/test_metrics.R

.args <- commandArgs(trailingOnly = FALSE)
.script <- sub("^--file=", "", .args[grep("^--file=", .args)])
.script_dir <- if (length(.script)) dirname(.script) else "benchmarks/study/common"
source(file.path(.script_dir, "metrics.R"))

fail_count <- 0L
check <- function(desc, got, want, tol = 1e-9) {
  ok <- if (is.infinite(want)) is.infinite(got) && sign(got) == sign(want)
        else isTRUE(all.equal(got, want, tolerance = tol))
  cat(sprintf("[%s] %-60s got=%s want=%s\n", if (ok) "PASS" else "FAIL",
              desc, format(got, digits = 12), format(want, digits = 12)))
  if (!ok) fail_count <<- fail_count + 1L
  invisible(ok)
}
check_true <- function(desc, got) {
  ok <- isTRUE(got)
  cat(sprintf("[%s] %-60s got=%s\n", if (ok) "PASS" else "FAIL", desc, got))
  if (!ok) fail_count <<- fail_count + 1L
  invisible(ok)
}

## ---------------------------------------------------------------------------
## Toy: 3x4 -- margin A (3 levels) x margin B (4 levels), n=6 observations.
## ---------------------------------------------------------------------------
groupA <- c("a1", "a1", "a2", "a2", "a3", "a3")
groupB <- c("b1", "b2", "b3", "b4", "b1", "b2")
TA <- c(a1 = 0.3, a2 = 0.3, a3 = 0.4)
TB <- c(b1 = 0.25, b2 = 0.25, b3 = 0.25, b4 = 0.25)
groups  <- list(A = groupA, B = groupB)
targets <- list(A = TA, B = TB)

## ---- Legacy reference impls (verbatim copies, cited; NOT modified in situ) -

# Verbatim from benchmarks/stepstone_all_methods.R:22-42 (fit_metrics), adapted
# only to take (w, groupA, groupB, TA, TB) instead of a data.frame + tgt list
# (same algorithm; toy has 2 margins named A/B instead of arbitrary df cols).
legacy_fit_metrics <- function(w, groupA, groupB, TA, TB) {
  W <- sum(w); n <- length(w)
  max_err <- 0; L1 <- 0; chi2 <- 0; marg_kl <- 0
  for (grp_tgt in list(list(groupA, TA), list(groupB, TB))) {
    grp <- grp_tgt[[1]]; tgt <- grp_tgt[[2]]
    lv  <- names(tgt)
    S   <- tapply(w, grp, sum, default = 0)[lv]; S[is.na(S)] <- 0
    tk  <- unname(tgt); Sr <- S / W
    max_err <- max(max_err, max(abs(Sr - tk)))
    L1      <- L1 + sum(abs(Sr - tk))
    exp_    <- tk * W
    chi2    <- chi2 + sum(ifelse(exp_ > 0, (S - exp_)^2 / exp_, 0))
    safe    <- tk > 0 & Sr > 0
    marg_kl <- marg_kl + sum(ifelse(safe, tk * log(tk / pmax(Sr, 1e-300)), 0))
  }
  pos <- w > 0
  weight_kl <- sum(w[pos] * log(w[pos])) / n   # BUG: hardcodes d_i=1
  list(max_err = max_err, L1 = L1, chi2 = chi2, marg_kl = marg_kl, weight_kl = weight_kl)
}

# Verbatim from benchmarks/allmethod_bench.R:30-50 (compute_metrics), adapted
# only to take (w, groupA, groupB, TA, TB) instead of a data.frame + tgt list.
legacy_allmethod_metrics <- function(w, groupA, groupB, TA, TB) {
  W <- sum(w); n <- length(w)
  max_err <- 0; mkl <- 0; wkl <- 0; chi2 <- 0
  for (grp_tgt in list(list(groupA, TA), list(groupB, TB))) {
    grp <- grp_tgt[[1]]; tgt <- grp_tgt[[2]]
    S <- tapply(w, grp, sum); S <- S / W
    T <- tgt[names(S)]
    T[is.na(T)] <- 0
    err <- abs(S - T)
    max_err <- max(max_err, max(err, na.rm = TRUE))
    pos <- T > 0 & S > 0
    mkl <- mkl + sum(T[pos] * log(T[pos] / S[pos]))
    chi2 <- chi2 + sum((S - T)^2 / pmax(T, 1e-12))
  }
  wm <- W / n
  pos <- w > 0
  wkl <- sum(w[pos] * log(w[pos] / wm)) / W    # BUG: hardcodes d_i=1 (via wm, not real d)
  c(max_err = max_err, marginal_kl = mkl, kl = wkl, chi2 = chi2)
}

## ---------------------------------------------------------------------------
## Case 0: baseline (uniform w, d_i=1) -- sanity + legacy agreement on marg_kl
## ---------------------------------------------------------------------------
cat("== Case 0: baseline uniform toy ==\n")
w0 <- c(1, 1, 1, 1, 1, 1)
m0 <- compute_metrics(w0, groups, targets, family = "kl")

check("marg_kl_mean (base)", m0$marg_kl_mean, 0.03430191557553895)
check("marg_kl_max (base)",  m0$marg_kl_max,  0.05889151782819174)
check("margin_linf (base)",  m0$margin_linf,  0.08333333333333334)
check("margin_l1 (base)",    m0$margin_l1,    0.4666666666666667)
check("ESS (base, uniform)", m0$ESS, 6)
check("DEFF_kish (base, uniform)", m0$DEFF_kish, 1)
check("weight_kl (base, d=1, uniform)", m0$weight_kl, 0)
check_true("g_weighted FALSE when d all 1", !m0$g_weighted)

# Legacy impls report Sum_k(KL_k) (not mean/max per margin); new impl's
# marg_kl_mean*2 == Sum_k(KL_k) for 2 margins -- both must equal the
# independently-derived total (klA_base + klB_base) on this non-divergent case.
total_kl_base <- 0.009712313322886162 + 0.05889151782819174
check("new impl total (marg_kl_mean*2) matches independent total (base)",
      m0$marg_kl_mean * 2, total_kl_base)
lg0 <- legacy_fit_metrics(w0, groupA, groupB, TA, TB)
check("legacy fit_metrics marg_kl matches independent total (base, non-divergent)",
      lg0$marg_kl, total_kl_base)
la0 <- legacy_allmethod_metrics(w0, groupA, groupB, TA, TB)
check("legacy allmethod_bench marginal_kl matches independent total (base, non-divergent)",
      unname(la0["marginal_kl"]), total_kl_base)

## ---------------------------------------------------------------------------
## Case (a): starved category -- assert marg_kl -> +Inf, divergent=TRUE, and
## that the legacy `t>0 & S>0` gate SILENTLY DROPS the divergent term (bug).
## ---------------------------------------------------------------------------
cat("\n== Case (a): starved category (b4 starved to w=0) ==\n")
w_starved <- c(1, 1, 1, 0, 1, 2)   # Sum(w)=6=n; b4's sole member (obs4) starved
m_starved <- compute_metrics(w_starved, groups, targets, family = "kl")

check("marg_kl_max is +Inf on starved category", m_starved$marg_kl_max, Inf)
check("marg_kl_mean is +Inf on starved category (propagates, not averaged away)",
      m_starved$marg_kl_mean, Inf)
check_true("divergent flag set TRUE", m_starved$marg_kl_divergent)
check_true("divergent_margins == 'B' (the starved margin)",
           identical(m_starved$marg_kl_divergent_margins, "B"))
ms_starved <- margin_stats(w_starved, groups, targets)
check_true("margin A stays finite (only B starved)", is.finite(ms_starved$per_margin$A$kl))
check("margin A KL value (starved case, still finite)", ms_starved$per_margin$A$kl, 0.05547042424760395)

# Confirm the LEGACY gate silently drops the divergent term (bug reproduction,
# not a golden failure of the new impl): legacy marg_kl must be FINITE here,
# proving DESIGN.md's "t>0 & S>0 gate flatters catastrophic failures" claim.
lg_starved <- legacy_fit_metrics(w_starved, groupA, groupB, TA, TB)
check_true("BUG REPRODUCED: legacy fit_metrics marg_kl is finite (silently drops divergent term)",
           is.finite(lg_starved$marg_kl))
la_starved <- legacy_allmethod_metrics(w_starved, groupA, groupB, TA, TB)
check_true("BUG REPRODUCED: legacy allmethod_bench marginal_kl is finite (silently drops divergent term)",
           is.finite(unname(la_starved["marginal_kl"])))

## ---------------------------------------------------------------------------
## Case (b): d_i != 1 -- assert g-weight DEFF != Kish DEFF, and that legacy
## weight_kl (hardcoded d_i=1) diverges from the new d-aware weight_kl (bug).
## ---------------------------------------------------------------------------
cat("\n== Case (b): non-uniform design weights d_i != 1 ==\n")
w_b <- c(1, 1, 1, 1, 2, 2)
d_b <- c(1, 1, 3, 3, 1, 1)
m_b <- compute_metrics(w_b, groups, targets, d = d_b, family = "kl")

check("ESS (case b)", m_b$ESS, 16 / 3)
check("DEFF_kish (case b, raw w)", m_b$DEFF_kish, 1.125)
check("DEFF_g (case b, g=w/d)", m_b$DEFF_g, 1.38)
check_true("g-weight DEFF != Kish DEFF (Gap D)", abs(m_b$DEFF_g - m_b$DEFF_kish) > 1e-6)
check_true("g_weighted TRUE when d_i != 1", m_b$g_weighted)
check("weight_kl (case b, real d_i)", m_b$weight_kl, 0.5753641449035616)

# Legacy impls hardcode d_i=1 -- their weight_kl on case-b's w must differ from
# the correct d-aware value (bug reproduction).
lg_b <- legacy_fit_metrics(w_b, groupA, groupB, TA, TB)
check("BUG REPRODUCED: legacy fit_metrics weight_kl (wrongly assumes d_i=1)",
      lg_b$weight_kl, 0.46209812037329684)
check_true("legacy weight_kl (d_i=1 assumed) != correct d-aware weight_kl",
           abs(lg_b$weight_kl - m_b$weight_kl) > 1e-3)

## bound violation on case-b weights, L=0, U=1.5
bv <- compute_metrics(w_b, groups, targets, d = d_b, bounds = list(L = 0, U = 1.5), family = "kl")
check("bound_viol_count", bv$bound_viol_count, 2)
check("bound_viol_max", bv$bound_viol_max, 0.5)
check("bound_viol_mean", bv$bound_viol_mean, 1 / 6)

## ---------------------------------------------------------------------------
## Case (c): minimax objective-value agreement (Blocker G) -- weight-vector
## correlation must NOT be computed for the minimax family.
## ---------------------------------------------------------------------------
cat("\n== Case (c): minimax objective-value agreement (Blocker G) ==\n")
w_mm     <- c(0.8, 1.1, 0.9, 1.3, 0.95, 0.95)
w_mm_ref <- c(0.7, 1.2, 1.0, 1.2, 0.95, 0.95)   # deliberately different vector, same L-inf optimum face
ag <- agreement(w_mm, w_mm_ref, family = "minimax",
                obj_val = 0.123456, obj_val_ref = 0.123457, obj_tol = 1e-5)
check_true("minimax agreement mode == 'objective_value' (not weight_vector)",
           identical(ag$mode, "objective_value"))
check("minimax abs_diff", ag$abs_diff, 1.000000000001e-06, tol = 1e-12)
check_true("minimax agree TRUE within obj_tol", ag$agree)
check_true("minimax agreement record has NO pearson/spearman fields (Blocker G)",
           is.null(ag$pearson) && is.null(ag$spearman))

ag_disagree <- agreement(w_mm, w_mm_ref, family = "minimax",
                          obj_val = 0.05, obj_val_ref = 0.20, obj_tol = 1e-6)
check_true("minimax agree FALSE when objective values differ beyond tol", !ag_disagree$agree)

# Strictly-convex family (kl) DOES get weight-vector stats (non-degenerate
# vectors -- w0 is uniform, sd=0, which would make cor() emit an NaN warning):
ag_kl <- agreement(w_b, c(1.1, 0.9, 1.0, 1.0, 2.05, 1.95), family = "kl")
check_true("kl-family agreement mode == 'weight_vector'", identical(ag_kl$mode, "weight_vector"))
check_true("kl-family agreement has pearson field", !is.null(ag_kl$pearson))

## ---------------------------------------------------------------------------
cat(sprintf("\n%s: %d assertion(s) failed.\n", if (fail_count == 0) "GOLDEN PASS" else "GOLDEN FAIL", fail_count))
quit(status = if (fail_count == 0) 0L else 1L, save = "no")
