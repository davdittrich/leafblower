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

fit_metrics <- function(w, df, tgt) {
  W <- sum(w); n <- length(w)
  max_err <- 0; L1 <- 0; chi2 <- 0; KL <- 0
  for (nm in names(tgt)) {
    lv <- names(tgt[[nm]])
    S  <- tapply(w, df[[nm]], sum, default=0)[lv]; S[is.na(S)] <- 0
    tk <- unname(tgt[[nm]]); Sr <- S/W
    max_err <- max(max_err, max(abs(Sr - tk)))
    L1   <- L1   + sum(abs(Sr - tk))
    exp_ <- tk * W
    chi2 <- chi2 + sum(ifelse(exp_>0, (S-exp_)^2/exp_, 0))
    safe <- tk>0 & Sr>0
    KL   <- KL   + sum(ifelse(safe, tk*log(tk/pmax(Sr,1e-300)), 0))
  }
  list(max_err=max_err, L1=L1, chi2=chi2, KL=KL,
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
  m    <- fit_metrics(w, df, tgt)
  iters  <- if (!is.null(res$iterations)) res$iterations else NA
  status <- if (!is.null(res$status)) res$status else 0L
  cat(sprintf("  wall=%6.1fs  iters=%4s  status=%d  max_err=%.4e  L1=%.4e  chi2=%.3e  KL=%.3e  DEFF=%.4f  ESS=%s  wmin=%.3f  wmax=%.3f\n",
    wall, ifelse(is.na(iters),"  —",iters), status,
    m$max_err, m$L1, m$chi2, m$KL, m$DEFF,
    format(round(m$ESS), big.mark=","), m$wmin, m$wmax))
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

cat("\n=== Pearson r vs iEPPA ===\n")
for (nm in c("raking","sinkhorn","grake","greg","cheby","autumn")) {
  rv <- get(paste0("r_", nm))
  cat(sprintf("  ieppa ↔ %-10s  r=%.6f\n", nm, cor(r_ieppa$w, rv$w)))
}
