#!/usr/bin/env Rscript
tryCatch(
  suppressPackageStartupMessages({
    library(leafblower)
    library(jsonlite)
  }),
  error = function(e) {
    cat("SKIP:", conditionMessage(e), "\n", file=stderr())
    quit(status=2, save="no")
  }
)
args <- commandArgs(trailingOnly=TRUE)
df  <- read.csv(args[1], stringsAsFactors=FALSE)
tgt <- fromJSON(args[2])
tgt <- lapply(tgt, function(x) { v <- unlist(x); v / sum(v) })
res <- harvest(df, tgt,
               method          = args[3],
               min_weight      = 0,
               max_weight      = 5,
               max_iterations  = as.integer(args[5]),
               convergence     = list(metric = "max_err",
                                      rule   = "improvement",
                                      tol    = 1e-3),
               verbose         = 0L)
w <- res[["weights"]]
write.csv(data.frame(weight = w), args[4], row.names = FALSE)
