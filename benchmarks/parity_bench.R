#!/usr/bin/env Rscript
# benchmarks/parity_bench.R
# R side of R-Python parity benchmark on stepstone fulldata.
# All methods, defaults (per-method natural metric), tol=1e-4.
# Saves per-method weights to feather for cross-language comparison.

suppressPackageStartupMessages({
  library(arrow); library(jsonlite); library(leafblower)
})
`%||%` <- function(a, b) if (!is.null(a)) a else b

DATA    <- "benchmarks/stepstone_fulldata_bench_data.parquet"
TARGETS <- "benchmarks/stepstone_fulldata_bench_targets.json"
WDIR    <- "benchmarks/results/parity/weights_R"
OUT_CSV <- "benchmarks/results/parity/parity_bench_R.csv"
MAX_ITER <- 3000L; MAX_WT <- 5; TOL <- 1e-4
NREPS   <- 3L
DEF     <- list(tol = TOL)   # defaults: per-method natural metric (harvest.R switch)

dir.create(WDIR, recursive = TRUE, showWarnings = FALSE)

cat("Loading data...\n")
df_raw <- as.data.frame(read_parquet(DATA))
tgt    <- lapply(fromJSON(TARGETS), function(x) { x <- unlist(x); x / sum(x) })
df     <- df_raw[, names(tgt), drop = FALSE]
for (v in names(tgt)) df[[v]] <- factor(df[[v]], levels = names(tgt[[v]]))
n <- nrow(df)
cat(sprintf("n=%s | margins=%d | max_iter=%d | tol=%.0e\n\n",
            format(n, big.mark = ","), length(tgt), MAX_ITER, TOL))

STATUS <- c("0" = "converged", "1" = "no_conv", "2" = "infeasible",
            "3" = "bad_arg",   "4" = "budget",  "5" = "stall")

compute_metrics <- function(w, df, tgt) {
  W <- sum(w)
  max_err <- 0; mkl <- 0; chi2 <- 0
  for (v in names(tgt)) {
    S <- tapply(w, df[[v]], sum); S <- S / W
    T <- tgt[[v]][names(S)]; T[is.na(T)] <- 0
    err <- abs(S - T)
    max_err <- max(max_err, max(err, na.rm = TRUE))
    pos <- T > 0 & S > 0
    mkl <- mkl + sum(T[pos] * log(T[pos] / S[pos]))
    chi2 <- chi2 + sum((S - T)^2 / pmax(T, 1e-12))
  }
  wm <- W / length(w)
  pos <- w > 0
  wkl <- sum(w[pos] * log(w[pos] / wm)) / W
  c(max_err = max_err, marginal_kl = mkl, kl = wkl, chi2 = chi2)
}

slug <- function(s) gsub("[^A-Za-z0-9]+", "_", s)

run <- function(label, conv, ...) {
  times <- numeric(NREPS)
  res_last <- NULL
  for (rep in seq_len(NREPS)) {
    t0 <- proc.time()["elapsed"]
    res <- suppressWarnings(tryCatch(
      harvest(df, tgt, min_weight = 0, max_weight = MAX_WT,
              max_iterations = MAX_ITER, convergence = conv, verbose = 0L,
              attach_weights = FALSE, ...),
      error = function(e) structure(list(), class = "error", message = e$message)
    ))
    times[rep] <- round((proc.time()["elapsed"] - t0) * 1000)
    res_last <- res
  }
  if (inherits(res_last, "error")) {
    cat(sprintf("  ERROR: %s\n", attr(res_last, "message")))
    return(NULL)
  }
  ri <- attr(res_last, "result")
  w  <- as.numeric(res_last)        # attach_weights=FALSE → numeric vector
  m  <- compute_metrics(w, df, tgt)
  st <- STATUS[as.character(ri$status)]; if (is.na(st)) st <- as.character(ri$status)

  # Save weights for cross-language comparison
  fpath <- file.path(WDIR, paste0(slug(label), ".feather"))
  write_feather(data.frame(w = w), fpath)

  data.frame(method = label,
             median_ms = median(times), min_ms = min(times), max_ms = max(times),
             status = st, iterations = ri$iterations, algorithm = ri$algorithm,
             max_err = m["max_err"], marginal_kl = m["marginal_kl"],
             kl = m["kl"], chi2 = m["chi2"],
             stringsAsFactors = FALSE)
}

cfg <- list(
  list("oris",                  DEF, method = "oris"),
  list("oris_accel",            DEF, method = "oris",      accelerate = TRUE),
  list("oris_greedy",           DEF, method = "oris",      scheduler = "greedy"),
  list("oris_greedy_accel",     DEF, method = "oris",      scheduler = "greedy", accelerate = TRUE),
  list("oris_soft",             DEF, method = "oris_soft"),
  list("oris_soft_accel",       DEF, method = "oris_soft", accelerate = TRUE),
  list("raking",                 DEF, method = "raking"),
  list("greenkhorn",             DEF, method = "greenkhorn"),
  list("greenkhorn_greedy",      DEF, method = "greenkhorn", scheduler = "greedy"),
  list("greenkhorn_accel",       DEF, method = "greenkhorn", accelerate = TRUE),
  list("sinkhorn",               DEF, method = "sinkhorn"),
  list("chebyshev",              DEF, method = "chebyshev"),
  list("greg",                   DEF, method = "greg"),
  list("logit",                  DEF, method = "logit"),
  list("newton_kl",              DEF, method = "newton_kl")
)

results <- vector("list", length(cfg))
for (i in seq_along(cfg)) {
  args <- cfg[[i]]; label <- args[[1]]; conv <- args[[2]]; args <- args[-(1:2)]
  cat(sprintf("[%2d/%2d] %-22s", i, length(cfg), label)); flush.console()
  results[[i]] <- do.call(run, c(list(label = label, conv = conv), args))
  if (!is.null(results[[i]])) {
    r <- results[[i]]
    cat(sprintf("  %7dms  %-10s  iters=%d  mkl=%.4g  max_err=%.4g\n",
                r$median_ms, r$status, r$iterations, r$marginal_kl, r$max_err))
  }
}

tab <- do.call(rbind, Filter(Negate(is.null), results))
write.csv(tab, OUT_CSV, row.names = FALSE)
cat(sprintf("\nSaved metrics: %s\nSaved weights: %s/*.feather\n", OUT_CSV, WDIR))
