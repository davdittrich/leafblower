library(leafblower)
set.seed(3); n <- 5000L
df  <- data.frame(v1=factor(sample(5, n, TRUE)))
tgt <- list(v1=setNames(c(0.4, 0.3, 0.15, 0.1, 0.05), as.character(1:5)))
r   <- harvest(df, tgt, method="oris",
               max_weight=1.8, min_weight=0, max_iterations=500,
               convergence=list(improvement=1e-4),
               attach_weights=FALSE)
saveRDS(list(df=df, tgt=tgt, max_weight=1.8, min_weight=0,
             max_iterations=500L,
             convergence=list(improvement=1e-4),
             weights=as.numeric(r),
             result=attr(r,"result")),
        "tests/testthat/fixtures/oris_pre_alm_ref.rds")
cat("Fixture written. n_weights:", length(as.numeric(r)),
    "max_error:", attr(r,"result")$max_error, "\n")
