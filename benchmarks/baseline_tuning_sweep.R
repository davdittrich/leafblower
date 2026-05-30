library(leafblower)
stopifnot(requireNamespace("arrow", quietly = TRUE))
data   <- arrow::read_parquet("benchmarks/stepstone_fulldata_bench_data.parquet")
target <- lapply(
  jsonlite::fromJSON("benchmarks/stepstone_fulldata_bench_targets.json"),
  function(x) { v <- unlist(x); v / sum(v) }
)

max_err_of <- function(w, data, target) {
  errs <- numeric(length(target))
  for (k in seq_along(target)) {
    tab <- tapply(w, data[[names(target)[k]]], sum) / sum(w)
    errs[k] <- max(abs(tab - target[[k]]))
  }
  max(errs)
}

grid <- expand.grid(
  bounds_mode    = c("cell", "unit"),
  max_iterations = c(3000L, 6000L, 10000L),
  stringsAsFactors = FALSE
)

results <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  g <- grid[i, , drop = FALSE]
  cat(sprintf("Running: bounds_mode=%s max_iterations=%d ...\n",
              g$bounds_mode, g$max_iterations))
  t0 <- Sys.time()
  w <- leafblower::harvest(
    data, target,
    max_weight     = 5,
    method         = "oris",
    bounds_mode    = g$bounds_mode,
    max_iterations = g$max_iterations,
    convergence    = list(absolute = 1e-10),
    attach_weights = FALSE
  )
  wall <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  data.frame(g, wall_s = wall, max_err = max_err_of(w, data, target))
}))

print(results, digits = 4)
saveRDS(results, "benchmarks/baseline_tuning_sweep.rds")

hit <- any(results$max_err <= 1.60e-3)
if (hit) {
  stop("FALSIFIED: existing knobs reach 1.60e-3 on stepstone-fulldata. ",
       "Overlays may be unnecessary. Halt and revisit plan.")
} else {
  cat("Confirmed: no existing-knob config reaches 1.60e-3. Overlays justified.\n")
}
