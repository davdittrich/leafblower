#!/usr/bin/env Rscript
# All-methods comparison on stepstone fulldata (1.58M rows, 9 margins).
# Runs each method sequentially. OMP_NUM_THREADS=1 set by caller.

suppressPackageStartupMessages({
  library(arrow); library(jsonlite); library(leafblower); library(autumn)
})

df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
df$uuid <- NULL
tgt <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
tgt <- lapply(tgt, function(t) { t <- unlist(t); t / sum(t) })
for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

cat(sprintf("n=%s  K=%d  margins: %s\n\n",
  format(nrow(df), big.mark=","), length(tgt),
  paste(names(tgt), collapse=", ")))

fit_metrics <- function(w, df, tgt, res = NULL) {
  W <- sum(w); n <- length(w)
  max_err <- 0; L1 <- 0; chi2 <- 0; marg_kl <- 0
  for (nm in names(tgt)) {
    lv <- names(tgt[[nm]])
    S  <- tapply(w, df[[nm]], sum, default=0)[lv]; S[is.na(S)] <- 0
    tk <- unname(tgt[[nm]]); Sr <- S/W
    max_err  <- max(max_err, max(abs(Sr - tk)))
    L1       <- L1 + sum(abs(Sr - tk))
    exp_     <- tk * W
    chi2     <- chi2 + sum(ifelse(exp_>0, (S-exp_)^2/exp_, 0))
    safe     <- tk>0 & Sr>0
    marg_kl  <- marg_kl + sum(ifelse(safe, tk*log(tk/pmax(Sr,1e-300)), 0))
  }
  weight_kl <- if (!is.null(res) && !is.null(res$convergence_used$solver_objective))
                 res$convergence_used$solver_objective else NA_real_
  list(max_err=max_err, L1=L1, chi2=chi2, marg_kl=marg_kl, weight_kl=weight_kl,
       DEFF=n*sum(w^2)/W^2, ESS=W^2/sum(w^2),
       wmin=min(w), wmed=median(w), wmax=max(w))
}

run <- function(method, ...) {
  cat(sprintf("%-12s ...", method)); flush.console()
  t0 <- proc.time()["elapsed"]
  r  <- suppressWarnings(
    if (method == "autumn")
      autumn::harvest(df, tgt, max_weight=5, min_weight=0, accelerate=TRUE, ...)
    else
      leafblower::harvest(df, tgt, method=method,
        max_weight=5, min_weight=0, ..., attach_weights=FALSE, verbose=0)
  )
  wall <- proc.time()["elapsed"] - t0
  w <- tryCatch(as.numeric(r),
    error = function(e) {
      if (!is.null(r$weights)) as.numeric(r$weights)
      else as.numeric(unlist(r)[seq_len(nrow(df))])
    })
  res  <- attr(r, "result")
  m    <- fit_metrics(w, df, tgt, res = res)
  iters  <- if (!is.null(res$iterations)) res$iterations else NA
  status <- if (!is.null(res$status)) res$status else 0L
  cat(sprintf("  wall=%6.1fs  iters=%4s  status=%d  max_err=%.4e  L1=%.4e  chi2=%.3e  marg_kl=%.3e  weight_kl=%s  DEFF=%.4f  ESS=%s  wmin=%.3f  wmax=%.3f\n",
    wall, ifelse(is.na(iters),"  —",iters), status,
    m$max_err, m$L1, m$chi2, m$marg_kl,
    ifelse(is.na(m$weight_kl), "  —", sprintf("%.3e", m$weight_kl)),
    m$DEFF, format(round(m$ESS), big.mark=","), m$wmin, m$wmax))
  invisible(list(w=w, wall=wall, m=m))
}

ITERS <- 5000L
cat("=== leafblower methods ===\n")
r_ieppa    <- run("ieppa",    max_iterations=ITERS)
r_raking   <- run("raking",   max_iterations=ITERS)
r_sinkhorn <- run("sinkhorn", max_iterations=ITERS)
r_grake    <- run("grake",    max_iterations=ITERS)
r_greg     <- run("greg",     max_iterations=ITERS)
r_cheby    <- run("chebyshev",max_iterations=ITERS)

cat("\n=== autumn (reference) ===\n")
r_autumn   <- run("autumn")

cat("\n=== ieppa overlay combinations ===\n")

# P-B: greedy margin scheduler (max-residual priority)
run("ieppa",
    scheduler = "greedy",
    max_iterations = ITERS)

# P-A: progressive bound tightening (3 levels, k=5)
run("ieppa",
    homotopy_levels = 3L,
    homotopy_start_factor = 5.0,
    homotopy_end_factor = 1.0,
    max_iterations = ITERS)

# P-A + Tang-eta: combined dynamic schedule
run("ieppa",
    homotopy_levels = 3L,
    homotopy_start_factor = 5.0,
    homotopy_end_factor = 1.0,
    eta_schedule = "tang_dynamic",
    eta_schedule_power = 0.5,
    max_iterations = ITERS)

# P-A + P-B + Tang-eta: all combined
run("ieppa",
    scheduler = "greedy",
    homotopy_levels = 3L,
    homotopy_start_factor = 5.0,
    homotopy_end_factor = 1.0,
    eta_schedule = "tang_dynamic",
    eta_schedule_power = 0.5,
    max_iterations = ITERS)

cat("\n=== Pearson r vs iEPPA ===\n")
for (nm in c("raking","sinkhorn","grake","greg","cheby","autumn")) {
  rv <- get(paste0("r_", nm))
  cat(sprintf("  ieppa ↔ %-10s  r=%.6f\n", nm, cor(r_ieppa$w, rv$w)))
}

cat("\n=== Python IPF implementations ===\n")
cat("Note: Python methods have NO max_weight/min_weight bounds.\n")
cat("ipfn: DataFrame mode on unique cells (K=9, n=1.58M).\n")

py_json <- tryCatch({
  py_out <- system(
    "OMP_NUM_THREADS=1 python3 benchmarks/python_ipf_benchmark.py 2>/dev/null",
    intern = TRUE)
  if (length(py_out) == 0L) stop("no output")
  jsonlite::fromJSON(paste(py_out, collapse = "\n"))
}, error = function(e) {
  cat(sprintf("Python benchmark failed: %s\n", conditionMessage(e)))
  NULL
})

if (!is.null(py_json)) {
  for (i in seq_len(nrow(py_json))) {
    r <- py_json[i, ]
    note <- if (!is.null(r$note) && !is.na(r$note))
      sprintf("  [%s]", r$note) else ""
    cat(sprintf(
      "%-26s  wall=%6.1fs  iters=%4s  status=%d  max_err=%.4e  marg_kl=%.3e  weight_kl=%s  DEFF=%.4f  ESS=%s  wmin=%.3f  wmax=%.3f%s\n",
      r$method, r$wall,
      ifelse(is.na(r$iters), "  —", as.character(r$iters)),
      r$status, r$max_err, r$marg_kl,
      ifelse(is.na(r$weight_kl) | is.nan(r$weight_kl), "       —",
             sprintf("%.3e", r$weight_kl)),
      r$DEFF,
      format(round(r$ESS), big.mark = ","),
      r$wmin, r$wmax, note))
  }
}
