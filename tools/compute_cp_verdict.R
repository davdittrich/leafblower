#!/usr/bin/env Rscript
# WU-4: Compute CP verdict from benchmarks/research/results/cp_summary.csv per spec Sec 1
# Decision Rule. Writes research/cp_verdict.txt; closes WU-5+WU-6 via bd if wu5_skip=true.

df <- read.csv("benchmarks/research/results/cp_summary.csv", stringsAsFactors = FALSE)
stopifnot("kk1204_K20" %in% df$fixture)
kk <- df[df$fixture == "kk1204_K20", ]
stopifnot(nrow(kk) == 3L)
kk_max_err <- median(kk$max_err, na.rm = TRUE)
kk_wall_s  <- median(kk$wall_s, na.rm = TRUE)
ss <- df[df$fixture == "stepstone_K9", ]
ss_parity <- ss$stepstone_parity_ratio[1]

beta_last <- median(kk$rate_exponent_last,    na.rm = TRUE)
r2_last   <- median(kk$rate_R2_last,          na.rm = TRUE)
n_last    <- median(kk$n_fit_points_last,     na.rm = TRUE)
beta_ergo <- median(kk$rate_exponent_ergodic, na.rm = TRUE)
r2_ergo   <- median(kk$rate_R2_ergodic,       na.rm = TRUE)
n_ergo    <- median(kk$n_fit_points_ergodic,  na.rm = TRUE)

rate_ok_last <- (!is.na(r2_last) && r2_last >= 0.9 && !is.na(n_last) && n_last >= 30 &&
                 !is.na(beta_last) && beta_last <= -0.8)
rate_ok_ergo <- (!is.na(r2_ergo) && r2_ergo >= 0.9 && !is.na(n_ergo) && n_ergo >= 30 &&
                 !is.na(beta_ergo) && beta_ergo <= -1.0)
rate_ok      <- rate_ok_last || rate_ok_ergo
ss_ok        <- (!is.na(ss_parity) && ss_parity <= 1.5)

verdict <- {
  if (is.na(kk_max_err) || kk_max_err >= 1e-2 || kk_wall_s > 60) {
    "FAIL"
  } else if (kk_max_err < 1e-3 && kk_wall_s <= 30 && rate_ok && ss_ok) {
    "PASS"
  } else if (kk_max_err < 1e-3 && kk_wall_s <= 30 && rate_ok && !ss_ok) {
    "PASS-kk1204-specialist"
  } else if (kk_max_err < 1e-2 && kk_wall_s <= 60) {
    "PARTIAL"
  } else {
    "FAIL"
  }
}

wu5_skip <- verdict %in% c("PASS", "PASS-kk1204-specialist")
con <- file("research/cp_verdict.txt", "w")
writeLines(c(
  verdict,
  paste0("wu5_skip=", tolower(as.character(wu5_skip))),
  sprintf("kk1204_max_err_median=%.6e", kk_max_err),
  sprintf("kk1204_wall_s_median=%.4f", kk_wall_s),
  sprintf("stepstone_parity_ratio=%.4f", ss_parity),
  sprintf("rate_exponent_last=%.4f", beta_last),
  sprintf("rate_R2_last=%.4f", r2_last),
  sprintf("n_fit_points_last=%d", as.integer(n_last)),
  sprintf("rate_exponent_ergodic=%.4f", beta_ergo),
  sprintf("rate_R2_ergodic=%.4f", r2_ergo),
  sprintf("n_fit_points_ergodic=%d", as.integer(n_ergo)),
  sprintf("rate_ok_last=%s", as.character(rate_ok_last)),
  sprintf("rate_ok_ergodic=%s", as.character(rate_ok_ergo)),
  sprintf("stepstone_ok=%s", as.character(ss_ok))
), con = con)
close(con)
cat("WU-4 verdict:", verdict, "wu5_skip:", wu5_skip, "\n")
