#!/usr/bin/env Rscript
# benchmarks/allmethod_bench.R
# Full-method comparison on stepstone fulldata (1.58M rows, 9 margins).
# Method-specific convergence metrics. All loss values computed post-hoc from weights.

suppressPackageStartupMessages({
  library(arrow); library(jsonlite); library(leafblower)
})

DATA    <- "benchmarks/stepstone_fulldata_bench_data.parquet"
TARGETS <- "benchmarks/stepstone_fulldata_bench_targets.json"
MAX_ITER <- 3000L; MAX_WT <- 5; TOL <- 1e-3

cat("Loading data...\n")
df_raw <- as.data.frame(read_parquet(DATA))
tgt    <- lapply(fromJSON(TARGETS), function(x){x<-unlist(x);x/sum(x)})
df     <- df_raw[,names(tgt),drop=FALSE]
for(v in names(tgt)) df[[v]] <- factor(df[[v]], levels=names(tgt[[v]]))
n <- nrow(df)
cat(sprintf("n=%s | margins=%d | max_iter=%d | tol=%.0e\n\n",
            format(n,big.mark=","), length(tgt), MAX_ITER, TOL))

STATUS <- c("0"="converged","1"="no_conv","2"="infeasible","3"="bad_arg","4"="budget","5"="stall")

# Post-hoc metric computation from weights
compute_metrics <- function(w, df, tgt) {
  W <- sum(w)
  max_err <- 0; mkl <- 0; wkl <- 0; chi2 <- 0
  w0 <- rep(1, length(w))   # uniform baseline for weight-KL
  for(v in names(tgt)) {
    S <- tapply(w, df[[v]], sum); S <- S / W
    T <- tgt[[v]][names(S)]
    T[is.na(T)] <- 0
    err <- abs(S - T)
    max_err <- max(max_err, max(err, na.rm=TRUE))
    # marginal KL: Σ T*log(T/S)  (skip zero T)
    pos <- T > 0 & S > 0
    mkl <- mkl + sum(T[pos] * log(T[pos] / S[pos]))
    chi2 <- chi2 + sum((S - T)^2 / pmax(T, 1e-12))
  }
  # weight-KL: Σ w_i * log(w_i / w_mean)
  wm <- W / n
  pos <- w > 0
  wkl <- sum(w[pos] * log(w[pos] / wm)) / W
  c(max_err=max_err, marginal_kl=mkl, kl=wkl, chi2=chi2)
}

run <- function(label, conv, ...) {
  t0 <- proc.time()["elapsed"]
  res <- suppressWarnings(tryCatch(
    harvest(df, tgt, min_weight=0, max_weight=MAX_WT,
            max_iterations=MAX_ITER, convergence=conv, verbose=0L, ...),
    error=function(e) structure(list(), class="error", message=e$message)
  ))
  elapsed <- round((proc.time()["elapsed"]-t0)*1000)
  if(inherits(res,"error")) {
    cat(sprintf("  ERROR: %s\n", attr(res,"message")))
    return(NULL)
  }
  ri  <- attr(res,"result")
  w   <- res[["weights"]]  # harvest returns data.frame with appended weight column
  m   <- compute_metrics(w, df, tgt)
  st  <- STATUS[as.character(ri$status)]
  if(is.na(st)) st <- as.character(ri$status)
  data.frame(method=label, loss_fn=conv$metric, time_ms=elapsed, status=st,
             iterations=ri$iterations, algorithm=ri$algorithm,
             max_err=signif(m["max_err"],4), marginal_kl=signif(m["marginal_kl"],4),
             kl=signif(m["kl"],4), chi2=signif(m["chi2"],4),
             stringsAsFactors=FALSE)
}

KL   <- list(metric="marginal_kl", rule="threshold", tol=TOL)
MERR <- list(metric="max_err",     rule="threshold", tol=TOL)
WKL  <- list(metric="kl",         rule="threshold", tol=TOL)
CHI  <- list(metric="chi2",       rule="threshold", tol=TOL)

cfg <- list(
  list("ieppa",                  KL,   method="ieppa"),
  list("ieppa + accel",          KL,   method="ieppa",      accelerate=TRUE),
  list("ieppa + greedy",         KL,   method="ieppa",      scheduler="greedy"),
  list("ieppa + greedy + accel", KL,   method="ieppa",      scheduler="greedy", accelerate=TRUE),
  list("ieppa_soft (auto cp)",   KL,   method="ieppa_soft"),               # auto capacity_penalty = M_cell/n
  list("ieppa_soft + accel",    KL,   method="ieppa_soft", accelerate=TRUE),
  list("raking",                 MERR, method="raking"),   # raking minimizes weight-KL; conv check = max_err
  list("greenkhorn",             MERR, method="greenkhorn"),
  list("greenkhorn + greedy",    MERR, method="greenkhorn", scheduler="greedy"),
  list("greenkhorn + accel",     MERR, method="greenkhorn", accelerate=TRUE),
  list("sinkhorn",               WKL,  method="sinkhorn"),
  list("chebyshev",              CHI,  method="chebyshev"),
  list("greg",                   CHI,  method="greg"),
  list("logit",                  MERR, method="logit"),
  list("newton_kl",              WKL,  method="newton_kl")
)

results <- vector("list", length(cfg))
for(i in seq_along(cfg)) {
  args <- cfg[[i]]; label <- args[[1]]; conv <- args[[2]]; args <- args[-(1:2)]
  cat(sprintf("[%2d/%2d] %-36s", i, length(cfg), label)); flush.console()
  results[[i]] <- do.call(run, c(list(label=label,conv=conv),args))
  if(!is.null(results[[i]])) {
    r <- results[[i]]
    cat(sprintf("  %7dms  %-10s  mkl=%.4f  max_err=%.4f  iters=%d\n",
                r$time_ms, r$status, r$marginal_kl, r$max_err, r$iterations))
  }
}

tab <- do.call(rbind, Filter(Negate(is.null), results))
tab_s <- tab[order(tab$time_ms),]

cat("\n=== Results (sorted by time_ms) ===\n")
print(tab_s[,c("method","loss_fn","time_ms","status","iterations",
               "marginal_kl","max_err","kl","chi2","algorithm")],
      row.names=FALSE, digits=4)

outfile <- "benchmarks/results/allmethod_bench_r.csv"
write.csv(tab, outfile, row.names=FALSE)
cat(sprintf("\nSaved: %s\n", outfile))
