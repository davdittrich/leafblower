# Prefer fulldata if available (main repo); fall back to standard bench data in worktree.
src_candidates <- c(
  "benchmarks/stepstone_fulldata_bench_data.parquet",
  "benchmarks/stepstone_bench_data.parquet"
)
src <- Filter(file.exists, src_candidates)[1]
stopifnot(!is.na(src))

tg_candidates <- c(
  "benchmarks/stepstone_fulldata_bench_targets.json",
  "benchmarks/stepstone_bench_targets.json"
)
tg_src <- Filter(file.exists, tg_candidates)[1]
stopifnot(!is.na(tg_src))

data <- arrow::read_parquet(src)
set.seed(42)
n_sample <- min(10000L, nrow(data))
small <- data[sample.int(nrow(data), n_sample), , drop = FALSE]
dir.create("tests/testthat/fixtures", showWarnings = FALSE, recursive = TRUE)
arrow::write_parquet(small, "tests/testthat/fixtures/stepstone_small.parquet")
# fromJSON returns nested lists; convert each margin to a named numeric vector.
# Strip literal "NA" keys (JSON artefact from jsonlite encoding R NA levels),
# remove categories absent from the sample, then renormalise to sum-to-1.
target_raw <- jsonlite::fromJSON(tg_src, simplifyVector = TRUE)
target_list <- lapply(target_raw, function(x) {
  v <- if (is.list(x)) unlist(x) else x
  # Drop literal "NA" string key
  v <- v[names(v) != "NA"]
  v
})

# Restrict to categories present in sample and renormalise each margin.
# Each top-level name in target_list corresponds to a column in small.
target <- lapply(names(target_list), function(nm) {
  tgt <- target_list[[nm]]
  col <- small[[nm]]
  if (is.null(col)) return(tgt / sum(tgt))  # column not in data; keep as-is
  sv <- as.character(unique(col[!is.na(col)]))
  keep <- names(tgt) %in% sv
  if (!all(keep)) tgt <- tgt[keep]
  tgt / sum(tgt)
})
names(target) <- names(target_list)
saveRDS(target, "tests/testthat/fixtures/stepstone_small_targets.rds")
cat("Fixtures written (source:", src, ").\n")
