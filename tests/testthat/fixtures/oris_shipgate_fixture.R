library(leafblower)
set.seed(20260531)
n <- 5000L
K <- 8L
# K=8 binary margins: margins 1-4 sampled at 0.85/0.15, targeted at 0.15/0.85
# margins 5-8 sampled at 0.15/0.85, targeted at 0.85/0.15
df <- data.frame(
  m1 = sample(c("a","b"), n, replace=TRUE, prob=c(0.85,0.15)),
  m2 = sample(c("a","b"), n, replace=TRUE, prob=c(0.85,0.15)),
  m3 = sample(c("a","b"), n, replace=TRUE, prob=c(0.85,0.15)),
  m4 = sample(c("a","b"), n, replace=TRUE, prob=c(0.85,0.15)),
  m5 = sample(c("a","b"), n, replace=TRUE, prob=c(0.15,0.85)),
  m6 = sample(c("a","b"), n, replace=TRUE, prob=c(0.15,0.85)),
  m7 = sample(c("a","b"), n, replace=TRUE, prob=c(0.15,0.85)),
  m8 = sample(c("a","b"), n, replace=TRUE, prob=c(0.15,0.85))
)
tgt <- list(
  m1=c(a=0.15, b=0.85), m2=c(a=0.15, b=0.85),
  m3=c(a=0.15, b=0.85), m4=c(a=0.15, b=0.85),
  m5=c(a=0.85, b=0.15), m6=c(a=0.85, b=0.15),
  m7=c(a=0.85, b=0.15), m8=c(a=0.85, b=0.15)
)
# Baseline 1: no SOR
res_no_sor <- harvest(df, tgt, method="oris", sor=NULL,
                      max_weight=1000, min_weight=0, max_iterations=2000)
iters_no_sor <- attr(res_no_sor, "result")$iterations
cat("no_sor iters:", iters_no_sor, "\n")

# Baseline 2: fixed omega=1.5 (mode_id=1: always use omega_max=1.5)
res_fixed <- harvest(df, tgt, method="oris", sor=list(auto=TRUE, omega_mode_id=1L),
                     max_weight=1000, min_weight=0, max_iterations=2000)
iters_fixed <- attr(res_fixed, "result")$iterations
cat("fixed_sor iters:", iters_fixed, "\n")

# Stepstone mw=5 (ORIS fixed) - use small n=10000 from parquet if available, else -1
stepstone_iters <- -1L
parquet_path <- "tests/testthat/fixtures/stepstone_small.parquet"
if (file.exists(parquet_path)) {
  library(arrow)
  ss <- as.data.frame(arrow::read_parquet(parquet_path))
  if (nrow(ss) > 10000) ss <- ss[1:10000, ]
  ss_tgt_file <- "tests/testthat/fixtures/stepstone_small_targets.rds"
  if (file.exists(ss_tgt_file)) {
    ss_tgt <- readRDS(ss_tgt_file)
    ss_res <- tryCatch(
      harvest(ss, ss_tgt, method="oris", sor=list(auto=FALSE),
              max_weight=5, min_weight=0, max_iterations=500),
      error = function(e) NULL
    )
    if (!is.null(ss_res))
      stepstone_iters <- attr(ss_res, "result")$iterations
  }
}
cat("stepstone_mw5_fixed iters:", stepstone_iters, "\n")
