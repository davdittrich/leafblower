suppressPackageStartupMessages({ library(leafblower); library(arrow); library(jsonlite); library(RhpcBLASctl) })
Sys.setenv(OMP_NUM_THREADS = "1")
blas_set_num_threads(1L)
dyn.load("research/leafblower_research.so")
source("benchmarks/research/utils.R")  # build_design_matrix, check_memory, make_*

stopifnot(file.exists("benchmarks/stepstone_bench_data.parquet"))
stopifnot(file.exists("benchmarks/stepstone_bench_targets.json"))

git_sha <- system("git rev-parse --short HEAD", intern = TRUE)
omp_threads <- Sys.getenv("OMP_NUM_THREADS")
blas_threads <- blas_get_num_procs()

fit_rate <- function(tr, col) {
  if (!col %in% names(tr) || sum(!is.na(tr[[col]])) < 30) return(list(beta=NA, r2=NA, n=0))
  iter_col <- if ("iter_outer" %in% names(tr)) tr$iter_outer else tr$iter
  keep <- iter_col > 5 & iter_col < 0.9 * max(iter_col) & is.finite(tr[[col]]) & tr[[col]] > 0
  if (sum(keep) < 30) return(list(beta=NA, r2=NA, n=sum(keep)))
  m <- lm(log(tr[[col]][keep]) ~ log(iter_col[keep]))
  list(beta = unname(coef(m)[2]), r2 = summary(m)$r.squared, n = sum(keep))
}

summary_rows <- list()
fixtures <- list(t1_small = make_t1_small(seed=1L),
                 stepstone_K9 = make_stepstone_K9(),
                 kk1204_K20 = make_kk1204_K20(seed=1L))

dir.create("benchmarks/research/results", showWarnings=FALSE, recursive=TRUE)

for (fx_name in names(fixtures)) {
  fx <- fixtures[[fx_name]]
  check_memory(fx$n, fx$K, max_iter = 5000L)
  reps <- if (fx_name == "kk1204_K20") 3L else 1L
  for (rep in seq_len(reps)) {
    t0 <- Sys.time()
    r <- .Call("ipm_solve_R", fx$A_csr, fx$b, fx$d, fx$lo, fx$hi, 5000L, TRUE, 1L)
    wall_s <- as.numeric(Sys.time() - t0, units = "secs")
    tr <- as.data.frame(r$trace)
    if (nrow(tr) > 0) {
      colnames(tr) <- c("iter_outer", "iter_inner", "time_ms", "mu",
                        "max_err", "kkt_resid", "n_projected_dims", "alpha_step")
    }
    trace_path <- sprintf("benchmarks/research/results/ipm_%s_trace_rep%d.csv", fx_name, rep)
    write.csv(tr, trace_path, row.names = FALSE)

    fit_last <- fit_rate(tr, "max_err")

    summary_rows[[length(summary_rows)+1L]] <- data.frame(
      solver = "ipm", fixture = fx_name, rep = rep,
      status_code = r$status_code, status_msg = r$status_msg,
      converged = (r$status_code == 0L),
      max_err = if (nrow(tr)) tail(tr$max_err, 1L) else NA,
      max_err_ergodic = NA,
      wall_s = wall_s, n_iter = r$iterations,
      final_primal_resid = NA,
      final_residual = if ("kkt_resid" %in% names(tr) && nrow(tr)) tail(tr$kkt_resid, 1L) else NA,
      n_projected_dims_max = if ("n_projected_dims" %in% names(tr) && nrow(tr)) max(tr$n_projected_dims, na.rm = TRUE) else NA,
      rate_exponent_last = fit_last$beta, rate_R2_last = fit_last$r2, n_fit_points_last = fit_last$n,
      rate_exponent_ergodic = NA, rate_R2_ergodic = NA, n_fit_points_ergodic = NA,
      stepstone_parity_ratio = if (fx_name == "stepstone_K9" && nrow(tr)) tail(tr$max_err, 1L) / 1.13e-4 else NA,
      seed = 1L, omp_threads = omp_threads, blas_threads = blas_threads, git_sha = git_sha,
      stringsAsFactors = FALSE
    )
  }
}
write.csv(do.call(rbind, summary_rows), "benchmarks/research/results/ipm_summary.csv", row.names = FALSE)
cat("WU-6 IPM bench complete\n")
