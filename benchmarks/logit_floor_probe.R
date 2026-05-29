## benchmarks/logit_floor_probe.R
## Measure iters/max_err across kDeffFloor variants.
## Called with: Rscript benchmarks/logit_floor_probe.R <label>
## Ticket: leafblower-xwqy

suppressPackageStartupMessages(library(leafblower))
label <- commandArgs(trailingOnly=TRUE)[1]
if (is.na(label)) label <- "unknown"

run <- function(name, seed, df_maker, tgt_maker, min_w, max_w) {
  set.seed(seed)
  df  <- df_maker()
  tgt <- tgt_maker(df)
  r <- tryCatch(suppressWarnings(harvest(df, tgt, method="logit",
    min_weight=min_w, max_weight=max_w,
    max_iterations=50L, convergence=list(absolute=1e-4))),
    error=function(e) NULL)
  if (is.null(r)) return(list(name=name, iters=NA, max_err=NA, status=-1L))
  res <- attr(r,"result")
  list(name=name, iters=res$iterations, max_err=res$max_error, status=res$status)
}

results <- list(
  ## C_benign: K=2 T7 clone — fast converging control
  run("C_benign_K2_T7", 1L,
    df_maker  = function() data.frame(
      g1=factor(rep(c("a","b"),c(200L,300L))),
      g2=factor(c(rep("x",150L),rep("y",350L)))),
    tgt_maker = function(df) list(g1=c(a=0.5,b=0.5),g2=c(x=0.4,y=0.6)),
    min_w=0, max_w=5),

  ## C_medium: K=3 balanced tight — medium pace
  run("C_medium_K3_tight", 10L,
    df_maker  = function() data.frame(
      v1=factor(sample(1:4,3000,TRUE)),
      v2=factor(sample(1:3,3000,TRUE)),
      v3=factor(sample(1:3,3000,TRUE))),
    tgt_maker = function(df) lapply(df, function(f) {
      p <- rep(1/nlevels(f),nlevels(f)); names(p) <- levels(f); p }),
    min_w=0.5, max_w=2.0),

  ## C_tight: K=4 shifted empirical — 25-iter converging tight config
  run("C_tight_K4_shifted", 42L,
    df_maker  = function() {
      n<-5000L; K<-4L; nj<-4L
      df<-as.data.frame(lapply(seq_len(K), function(k)
        factor(sample(seq_len(nj),n,replace=TRUE))))
      names(df)<-paste0("v",seq_len(K)); df },
    tgt_maker = function(df) lapply(df, function(f) {
      p<-prop.table(table(f))
      p2<-p+c(.10,-.05,-.03,-.02)[seq_along(p)]; p2/sum(p2) }),
    min_w=0.5, max_w=2.0)
)

cat(sprintf("\n=== floor variant: %s ===\n", label))
out <- do.call(rbind, lapply(results, function(r) {
  cat(sprintf("  %-25s iters=%3d  max_err=%.3e  status=%d\n",
              r$name, r$iters, r$max_err, r$status))
  data.frame(config=r$name, floor=label,
             iters=r$iters, max_err=r$max_err, status=r$status,
             stringsAsFactors=FALSE)
}))

dir.create("benchmarks/results", showWarnings=FALSE, recursive=TRUE)
csv <- "benchmarks/results/logit_floor_probe.csv"
write.table(out, csv, sep=",", col.names=!file.exists(csv),
            append=file.exists(csv), row.names=FALSE, quote=TRUE)
invisible(out)
