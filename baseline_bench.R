suppressPackageStartupMessages(library(arrow))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(leafblower))

df <- arrow::read_parquet('benchmarks/stepstone_fulldata_bench_data.parquet')
df$uuid <- NULL

tgt <- lapply(
  jsonlite::fromJSON('benchmarks/stepstone_fulldata_bench_targets.json'),
  function(t) { t <- unlist(t); t/sum(t) }
)

for(nm in names(tgt)) df[[nm]] <- factor(df[[nm]])

r <- leafblower::harvest(df, tgt, method='raking', accelerate=TRUE,
      max_weight=5, max_iterations=5000L, attach_weights=FALSE, verbose=0)

w <- as.numeric(r)
W <- sum(w)
n <- length(w)

cat('baseline: max_err=', attr(r,'result')$max_error,
    ' DEFF=', n*sum(w^2)/W^2, '\n')
