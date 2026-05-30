#!/usr/bin/env Rscript
# benchmarks/allmethod_bench.R
# Full-method comparison on stepstone fulldata (1.58M rows, 9 margins).
# Method-specific convergence metrics. All loss values computed post-hoc from weights.

suppressPackageStartupMessages({
  library(arrow); library(jsonlite); library(leafblower)
})
`%||%` <- function(a, b) if (!is.null(a)) a else b

DATA    <- "benchmarks/stepstone_fulldata_bench_data.parquet"
TARGETS <- "benchmarks/stepstone_fulldata_bench_targets.json"
MAX_ITER <- 3000L; MAX_WT <- 5; TOL <- 1e-4
NREPS   <- 3L        # repetitions per method for median timing
# Default convergence per method (tol only; metric set by harvest.R per-method defaults)
DEF <- list(tol=TOL)  # uses each method's natural metric (harvest.R switch)

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
  times <- numeric(NREPS)
  res_last <- NULL
  for (rep in seq_len(NREPS)) {
    t0 <- proc.time()["elapsed"]
    res <- suppressWarnings(tryCatch(
      harvest(df, tgt, min_weight=0, max_weight=MAX_WT,
              max_iterations=MAX_ITER, convergence=conv, verbose=0L,
              attach_weights=FALSE, ...),
      error=function(e) structure(list(), class="error", message=e$message)
    ))
    times[rep] <- round((proc.time()["elapsed"]-t0)*1000)
    res_last <- res
  }
  if(inherits(res_last,"error")) {
    cat(sprintf("  ERROR: %s\n", attr(res_last,"message")))
    return(NULL)
  }
  ri  <- attr(res_last,"result")
  w   <- res_last  # attach_weights=FALSE: harvest returns numeric weights vector directly
  m   <- compute_metrics(w, df, tgt)
  st  <- STATUS[as.character(ri$status)]
  if(is.na(st)) st <- as.character(ri$status)
  loss_fn <- conv$metric %||% ri$convergence_used$metric %||% "default"
  data.frame(method=label, loss_fn=loss_fn,
             median_ms=median(times), min_ms=min(times), max_ms=max(times),
             status=st, iterations=ri$iterations, algorithm=ri$algorithm,
             max_err=signif(m["max_err"],4), marginal_kl=signif(m["marginal_kl"],4),
             kl=signif(m["kl"],4), chi2=signif(m["chi2"],4),
             stringsAsFactors=FALSE)
}

# Per-method natural default: metric set by harvest.R switch (1ab1165).
# Only tol is set here; each method uses its solver-appropriate objective.
cfg <- list(
  list("oris",                  DEF, method="oris"),
  list("oris + accel",          DEF, method="oris",      accelerate=TRUE),
  list("oris + greedy",         DEF, method="oris",      scheduler="greedy"),
  list("oris + greedy + accel", DEF, method="oris",      scheduler="greedy", accelerate=TRUE),
  list("oris_soft (auto cp)",   DEF, method="oris_soft"),
  list("oris_soft + accel",     DEF, method="oris_soft", accelerate=TRUE),
  list("raking",                 DEF, method="raking"),
  list("greenkhorn",             DEF, method="greenkhorn"),
  list("greenkhorn + greedy",    DEF, method="greenkhorn", scheduler="greedy"),
  list("greenkhorn + accel",     DEF, method="greenkhorn", accelerate=TRUE),
  list("sinkhorn",               DEF, method="sinkhorn"),
  list("chebyshev",              DEF, method="chebyshev"),
  list("greg",                   DEF, method="greg"),
  list("logit",                  DEF, method="logit"),
  list("newton_kl",              DEF, method="newton_kl")
)

results <- vector("list", length(cfg))
for(i in seq_along(cfg)) {
  args <- cfg[[i]]; label <- args[[1]]; conv <- args[[2]]; args <- args[-(1:2)]
  cat(sprintf("[%2d/%2d] %-36s", i, length(cfg), label)); flush.console()
  results[[i]] <- do.call(run, c(list(label=label,conv=conv),args))
  if(!is.null(results[[i]])) {
    r <- results[[i]]
    cat(sprintf("  %7dms [%d-%d]  %-10s  mkl=%.4f  max_err=%.4f  iters=%d\n",
                r$median_ms, r$min_ms, r$max_ms, r$status, r$marginal_kl, r$max_err, r$iterations))
  }
}

tab <- do.call(rbind, Filter(Negate(is.null), results))
tab_s <- tab[order(tab$median_ms),]

cat("\n=== Results (sorted by median_ms) ===\n")
print(tab_s[,c("method","loss_fn","median_ms","min_ms","max_ms","status","iterations",
               "marginal_kl","max_err","kl","chi2","algorithm")],
      row.names=FALSE, digits=4)

outfile <- "benchmarks/results/allmethod_bench_r.csv"
write.csv(tab, outfile, row.names=FALSE)
cat(sprintf("\nSaved: %s\n", outfile))
