#!/usr/bin/env Rscript
# Lean all-methods wall-time timer on stepstone fulldata for speed-regression checks.
# Robust to a method that INFEAS-halts (tryCatch). Single-thread BLAS set by caller.
suppressPackageStartupMessages({library(arrow); library(jsonlite); library(leafblower)})
df  <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet"); df$uuid <- NULL
tgt <- jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json")
tgt <- lapply(tgt, function(t){t <- unlist(t); t/sum(t)})
for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

time1 <- function(label, method, ...) {
  # warm one GC out of the timed region
  invisible(gc())
  t0 <- proc.time()["elapsed"]
  out <- tryCatch({
    r <- suppressWarnings(leafblower::harvest(df, tgt, method = method, max_weight = 5,
           min_weight = 0, max_iterations = 5000L, attach_weights = FALSE, verbose = 0, ...))
    res <- attr(r, "result")
    sprintf("wall=%6.2fs status=%d iters=%s", proc.time()["elapsed"] - t0,
            res$status, ifelse(is.null(res$iterations), "NA", res$iterations))
  }, error = function(e) sprintf("wall=%6.2fs ERROR: %s", proc.time()["elapsed"] - t0,
                                 sub("\n.*", "", conditionMessage(e))))
  cat(sprintf("%-18s %s\n", label, out)); flush.console()
}

cfgs <- list(
  c("oris","oris"), c("oris_soft","oris_soft"),
  c("raking","raking"), c("sinkhorn","sinkhorn"),
  c("greenkhorn","greenkhorn"), c("logit","logit"),
  c("greg","greg"), c("chebyshev","chebyshev"))
for (cf in cfgs) time1(cf[1], cf[2])
# accel variants that exercise SRAA in the changed greenkhorn/oris TUs
time1("oris+greedy",      "oris",       scheduler = "greedy")
time1("raking+sraa",      "raking",     accelerate = TRUE)
time1("greenkhorn+sraa",  "greenkhorn", accelerate = TRUE)
