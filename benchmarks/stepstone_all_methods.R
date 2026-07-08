#!/usr/bin/env Rscript
# All-methods comparison on stepstone fulldata (1.58M rows, 9 margins).
# Runs each method sequentially. OMP_NUM_THREADS=1 set by caller.

suppressPackageStartupMessages({
  library(arrow); library(jsonlite); library(leafblower)
  if (requireNamespace("autumn", quietly = TRUE)) library(autumn)
})

df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
df$uuid <- NULL
tgt <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
tgt <- lapply(tgt, function(t) { t <- unlist(t); t / sum(t) })
for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

cat(sprintf("n=%s  K=%d  margins: %s\n\n",
  format(nrow(df), big.mark=","), length(tgt),
  paste(names(tgt), collapse=", ")))

# All metrics computed from final weights — no reliance on solver diagnostics.
# weight_kl = Σ w[i]*log(w[i]/d[i]) / n  (d[i]=1: uniform normalized start).
fit_metrics <- function(w, df, tgt) {
  W <- sum(w); n <- length(w)
  max_err <- 0; L1 <- 0; chi2 <- 0; marg_kl <- 0
  for (nm in names(tgt)) {
    lv  <- names(tgt[[nm]])
    S   <- tapply(w, df[[nm]], sum, default = 0)[lv]; S[is.na(S)] <- 0
    tk  <- unname(tgt[[nm]]); Sr <- S / W
    max_err <- max(max_err, max(abs(Sr - tk)))
    L1      <- L1 + sum(abs(Sr - tk))
    exp_    <- tk * W
    chi2    <- chi2 + sum(ifelse(exp_ > 0, (S - exp_)^2 / exp_, 0))
    safe    <- tk > 0 & Sr > 0
    marg_kl <- marg_kl + sum(ifelse(safe, tk * log(tk / pmax(Sr, 1e-300)), 0))
  }
  pos       <- w > 0
  weight_kl <- sum(w[pos] * log(w[pos])) / n   # d_i=1 uniform start → log(w/d)=log(w)
  list(max_err = max_err, L1 = L1, chi2 = chi2, marg_kl = marg_kl,
       weight_kl = weight_kl,
       DEFF = n * sum(w^2) / W^2, ESS = W^2 / sum(w^2),
       wmin = min(w), wmed = median(w), wmax = max(w))
}

FMT <- "  wall=%7.1fs  iters=%5s  status=%d  max_err=%.4e  L1=%.4e  chi2=%.3e  marg_kl=%.3e  weight_kl=%.3e  DEFF=%.4f  ESS=%s  wmin=%.3f  wmed=%.3f  wmax=%.3f\n"

run <- function(method, label = method, ...) {
  cat(sprintf("%-22s ...", label)); flush.console()
  t0 <- proc.time()["elapsed"]
  # Guard: a method that stop()s (e.g. logit INFEAS) must not abort the whole
  # benchmark. Catch, print an ERROR line, return a NULL-weight sentinel, continue.
  r  <- tryCatch(suppressWarnings(
    if (method == "autumn")
      autumn::harvest(df, tgt, max_weight = 5, min_weight = 0, accelerate = TRUE, ...)
    else
      leafblower::harvest(df, tgt, method = method,
        max_weight = 5, min_weight = 0, ..., attach_weights = FALSE, verbose = 0)
  ), error = function(e) e)
  wall <- proc.time()["elapsed"] - t0
  if (inherits(r, "error")) {
    cat(sprintf("  wall=%7.1fs  ERROR: %s\n", wall,
                sub("\n.*", "", conditionMessage(r))))
    return(invisible(list(w = NULL, wall = wall, m = NULL)))
  }
  w <- tryCatch(as.numeric(r),
    error = function(e) {
      if (!is.null(r$weights)) as.numeric(r$weights)
      else as.numeric(unlist(r)[seq_len(nrow(df))])
    })
  res    <- attr(r, "result")
  m      <- fit_metrics(w, df, tgt)
  iters  <- if (!is.null(res$iterations)) res$iterations else NA
  status <- if (!is.null(res$status)) res$status else 0L
  cat(sprintf(FMT, wall, ifelse(is.na(iters), "   —", iters), status,
    m$max_err, m$L1, m$chi2, m$marg_kl, m$weight_kl,
    m$DEFF, format(round(m$ESS), big.mark = ","),
    m$wmin, m$wmed, m$wmax))
  invisible(list(w = w, wall = wall, m = m))
}

ITERS <- 5000L
cat("=== leafblower methods ===\n")
r_oris        <- run("oris",       max_iterations = ITERS)
r_oris_greedy <- run("oris",       label = "oris+greedy",     max_iterations = ITERS, scheduler = "greedy")
r_oris_soft   <- run("oris_soft",  max_iterations = ITERS)
r_raking       <- run("raking",      max_iterations = ITERS)
r_raking_sraa  <- run("raking",      label = "raking+sraa",      max_iterations = ITERS, accelerate = TRUE)
r_raking_greedy <- run("raking",     label = "raking+greedy",    max_iterations = ITERS, scheduler = "greedy")
r_sinkhorn     <- run("sinkhorn",    max_iterations = ITERS)
r_greenkhorn   <- run("greenkhorn",  max_iterations = ITERS)
r_grk_sraa     <- run("greenkhorn",  label = "greenkhorn+sraa",  max_iterations = ITERS, accelerate = TRUE)
r_logit        <- run("logit",       max_iterations = ITERS)
# lbfgsb removed — strictly dominated by newton_kl
r_greg         <- run("greg",        max_iterations = ITERS)
r_cheby        <- run("chebyshev",   max_iterations = ITERS)

cat("\n=== autumn (reference — cached) ===\n")
{
  ac       <- readRDS("tests/testthat/fixtures/stepstone_reference_autumn_only.rds")
  m        <- fit_metrics(ac$autumn_weights, df, tgt)
  cat(sprintf("%-22s", "autumn(cached)")); cat(" ...")
  cat(sprintf(FMT, ac$autumn_wall_s, "   —", 0L,
    m$max_err, m$L1, m$chi2, m$marg_kl, m$weight_kl,
    m$DEFF, format(round(m$ESS), big.mark = ","),
    m$wmin, m$wmed, m$wmax))
  r_autumn <- list(w = ac$autumn_weights, wall = ac$autumn_wall_s, m = m)
}

cat("\n=== Pearson r vs ORIS ===\n")
pairs <- list(
  oris_greedy  = "oris+greedy",
  oris_soft    = "oris_soft",
  raking        = "raking",
  raking_sraa   = "raking+sraa",
  raking_greedy = "raking+greedy",
  sinkhorn      = "sinkhorn",
  greenkhorn    = "greenkhorn",
  grk_sraa      = "greenkhorn+sraa",
  logit         = "logit",
  # lbfgsb removed — strictly dominated by newton_kl
  greg          = "greg",
  cheby         = "chebyshev",
  autumn        = "autumn"
)
if (is.null(r_oris$w)) {
  cat("  (reference oris run failed — Pearson-r section skipped)\n")
} else for (nm in names(pairs)) {
  rv <- get(paste0("r_", nm))
  if (is.null(rv$w) || length(rv$w) != length(r_oris$w)) {
    cat(sprintf("  oris ↔ %-22s  r=      — (failed)\n", pairs[[nm]])); next
  }
  cat(sprintf("  oris ↔ %-22s  r=%.6f\n", pairs[[nm]], cor(r_oris$w, rv$w)))
}

cat("\n=== Python IPF implementations ===\n")
cat("Note: Python methods have NO max_weight/min_weight bounds.\n")

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
    r    <- py_json[i, ]
    note <- if (!is.null(r$note) && !is.na(r$note)) sprintf("  [%s]", r$note) else ""
    cat(sprintf(
      "%-22s ...  wall=%7.1fs  iters=%5s  status=%d  max_err=%.4e  marg_kl=%.3e  weight_kl=%s  DEFF=%.4f  ESS=%s  wmin=%.3f  wmax=%.3f%s\n",
      r$method, r$wall,
      ifelse(is.na(r$iters), "   —", as.character(r$iters)),
      r$status, r$max_err, r$marg_kl,
      ifelse(is.na(r$weight_kl) | is.nan(r$weight_kl), "       —",
             sprintf("%.3e", r$weight_kl)),
      r$DEFF, format(round(r$ESS), big.mark = ","),
      r$wmin, r$wmax, note))
  }
}
