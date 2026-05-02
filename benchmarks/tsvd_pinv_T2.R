#!/usr/bin/env Rscript
# WL-2: T2 stepstone pinv-only measurement gate
Sys.setenv(OMP_NUM_THREADS = "1")
suppressPackageStartupMessages({ library(leafblower); library(arrow); library(jsonlite); library(bench) })

data_path <- {
  if (file.exists("benchmarks/stepstone_bench_data.parquet"))
    "benchmarks/stepstone_bench_data.parquet"
  else "../../benchmarks/stepstone_bench_data.parquet"
}
tgt_path <- {
  if (file.exists("benchmarks/stepstone_bench_targets.json"))
    "benchmarks/stepstone_bench_targets.json"
  else "../../benchmarks/stepstone_bench_targets.json"
}

df  <- read_parquet(data_path); df$uuid <- NULL
tgt <- lapply(jsonlite::fromJSON(tgt_path), function(t) { t <- unlist(t); t/sum(t) })
for (nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

t0 <- Sys.time()
r <- suppressWarnings(harvest(df, tgt, method="newton_kl",
  max_weight=5, min_weight=0, max_iterations=200L,
  attach_weights=FALSE, verbose=0L))
wall_ms <- as.numeric(Sys.time()-t0, units="secs") * 1000
res <- attr(r, "result")

result_row <- data.frame(
  gap = res$max_error,
  n_projected_dims = ifelse(is.null(res$n_projected_dims), NA_integer_, res$n_projected_dims),
  n_iter = res$iterations,
  wall_time_ms = round(wall_ms, 2)
)

dir.create("benchmarks/results", showWarnings=FALSE, recursive=TRUE)
write.csv(result_row, "benchmarks/results/tsvd_pinv_T2.csv", row.names=FALSE)

cat(sprintf("[wl-2] T2 stepstone K=9: gap=%.3e, n_proj=%d, n_iter=%d, wall=%.0fms\n",
  result_row$gap, result_row$n_projected_dims, result_row$n_iter, result_row$wall_time_ms))

if (result_row$gap < 1e-4) {
  cat("[wl-2] PASS — pinv alone clears T2 gate. CG retention TBD by WL-6.\n")
} else {
  cat(sprintf("[wl-2] MEASUREMENT — pinv insufficient (gap=%.3e ≥ 1e-4). CG required by WL-3. PROCEED.\n",
    result_row$gap))
}
