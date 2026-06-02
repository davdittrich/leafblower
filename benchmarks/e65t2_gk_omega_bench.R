# e65t.2 benchmark: greenkhorn gk_omega sweep
# Tests omega in {1.0, 1.2, 1.4} x 4 configs x 3 seeds
# Verdict: GO if omega>1 reduces iterations >=10% with zero divergence
library(leafblower)

set.seed(NULL)  # use per-config seeds explicitly

OMEGAS  <- c(1.0, 1.2, 1.4)
SEEDS   <- c(42L, 123L, 999L)

# 4 configs: K=3 small, K=6 medium, K=9 bounded, K=6 tight-bounds
configs <- list(
  list(name="K3_small",   n=200,  K=3, max_weight=3.0, min_weight=0.3),
  list(name="K6_medium",  n=500,  K=6, max_weight=3.0, min_weight=0.3),
  list(name="K9_bounded", n=800,  K=9, max_weight=2.0, min_weight=0.5),
  list(name="K6_tight",   n=400,  K=6, max_weight=1.5, min_weight=0.5)
)

run_one <- function(cfg, omega, seed) {
  set.seed(seed)
  n <- cfg$n; K <- cfg$K
  col_names <- paste0("x", seq_len(K))
  df <- data.frame(lapply(seq_len(K), function(k) {
    sample(LETTERS[1:4], n, replace=TRUE)
  }))
  names(df) <- col_names
  targets <- setNames(lapply(seq_len(K), function(k) {
    levs <- sort(unique(df[[col_names[k]]]))
    p    <- runif(length(levs)); p <- p / sum(p)
    setNames(p, levs)
  }), col_names)
  tryCatch({
    res <- harvest(df, targets,
                   method      = "greenkhorn",
                   gk_omega    = omega,
                   max_weight  = cfg$max_weight,
                   min_weight  = cfg$min_weight,
                   convergence = list(metric = "kl"),
                   accelerate  = FALSE,
                   max_iterations = 2000L)
    r <- attr(res, "result")
    list(iters=r$iterations, kl=r$kl, status=r$status, diverged=FALSE)
  }, error=function(e) {
    message("ERROR [", cfg$name, " omega=", omega, " seed=", seed, "]: ", conditionMessage(e))
    list(iters=NA_integer_, kl=NA_real_, status=-1L, diverged=TRUE)
  })
}

results <- list()
for (cfg in configs) {
  for (omega in OMEGAS) {
    for (seed in SEEDS) {
      key <- sprintf("%s_w%.1f_s%d", cfg$name, omega, seed)
      cat(sprintf("Running %s ...\n", key))
      r <- run_one(cfg, omega, seed)
      results[[key]] <- c(list(config=cfg$name, omega=omega, seed=seed), r)
    }
  }
}

# Summarize
res_df <- do.call(rbind, lapply(names(results), function(k) {
  r <- results[[k]]
  data.frame(config   = r$config,
             omega    = r$omega,
             seed     = r$seed,
             iters    = if (length(r$iters)   == 1L) r$iters   else NA_integer_,
             kl       = if (length(r$kl)      == 1L) r$kl      else NA_real_,
             status   = if (length(r$status)  == 1L) r$status  else NA_integer_,
             diverged = r$diverged,
             stringsAsFactors = FALSE)
}))

cat("\n=== RESULTS ===\n")
print(res_df)

# Verdict
diverge_count <- sum(res_df$diverged, na.rm=TRUE)
baseline <- subset(res_df, omega == 1.0 & !diverged)
omega12   <- subset(res_df, omega == 1.2 & !diverged)
omega14   <- subset(res_df, omega == 1.4 & !diverged)

# Median iteration improvement per omega (matched by config+seed)
med_imp_12 <- median((baseline$iters - omega12$iters) / baseline$iters * 100, na.rm=TRUE)
med_imp_14 <- median((baseline$iters - omega14$iters) / baseline$iters * 100, na.rm=TRUE)

cat(sprintf("\nDivergence count (omega>1): %d\n", diverge_count))
cat(sprintf("Median iter improvement omega=1.2: %.1f%%\n", med_imp_12))
cat(sprintf("Median iter improvement omega=1.4: %.1f%%\n", med_imp_14))

# GO: >=10% improvement at either omega level, zero divergence
verdict <- if (diverge_count == 0 && (med_imp_12 >= 10 || med_imp_14 >= 10)) "GO" else "NO-GO"
cat(sprintf("\nVERDICT: %s\n", verdict))
