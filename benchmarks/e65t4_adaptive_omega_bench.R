# e65t.4.5: Adaptive omega ship-gate benchmark
# 4 arms x 4 configs x 3 seeds; metric="kl"; accelerate=FALSE; max_iterations=2000L
library(leafblower)

set.seed(0)

configs <- list(
  list(name="unbounded_K3", n=200,  K=3, max_weight=3.0, min_weight=0.3),
  list(name="unbounded_K6", n=500,  K=6, max_weight=3.0, min_weight=0.3),
  list(name="bounded_K6",   n=500,  K=6, max_weight=2.0, min_weight=0.4),
  list(name="bounded_K9",   n=800,  K=9, max_weight=1.5, min_weight=0.5)
)
seeds <- c(42L, 123L, 999L)
arms <- list(
  list(name="baseline",  sk_omega=1.0, sk_omega_mode_id=0L),
  list(name="fixed_1.2", sk_omega=1.2, sk_omega_mode_id=0L),
  list(name="fixed_1.4", sk_omega=1.4, sk_omega_mode_id=0L),
  list(name="adaptive",  sk_omega=1.0, sk_omega_mode_id=1L)
)

# Helper: generate a random K-dimensional dataset with crossed factor
make_data <- function(n, K, seed) {
  set.seed(seed)
  df <- as.data.frame(
    lapply(seq_len(K), function(k) {
      nl <- 2L + (k %% 3L)  # 2-4 levels per margin
      sample(paste0("x", k, "_", seq_len(nl)), n, replace=TRUE)
    })
  )
  names(df) <- paste0("v", seq_len(K))
  df
}

make_targets <- function(df, K) {
  lapply(seq_len(K), function(k) {
    v <- df[[paste0("v", k)]]
    lvs <- sort(unique(v))
    p <- runif(length(lvs)); p <- p / sum(p)
    setNames(p, lvs)
  }) |> setNames(paste0("v", seq_len(K)))
}

results <- list()
for (cfg in configs) {
  for (seed in seeds) {
    df <- make_data(cfg$n, cfg$K, seed)
    targets <- make_targets(df, cfg$K)
    for (arm in arms) {
      raw <- tryCatch({
        harvest(
          df, targets,
          method          = "sinkhorn",
          max_iterations  = 2000L,
          max_weight      = cfg$max_weight,
          min_weight      = cfg$min_weight,
          sk_omega        = arm$sk_omega,
          sk_omega_mode_id = arm$sk_omega_mode_id,
          convergence     = list(metric = "kl"),
          accelerate      = FALSE
        )
      }, error = function(e) list(error=conditionMessage(e)))
      # harvest() returns a data.frame with attr "result" when attach_weights=TRUE
      res <- if (!is.null(attr(raw, "result"))) attr(raw, "result") else raw
      results[[length(results)+1]] <- list(
        config   = cfg$name,
        seed     = seed,
        arm      = arm$name,
        iters    = if(is.null(res$iterations)) NA_integer_ else as.integer(res$iterations),
        kl       = if(is.null(res$kl)) NA_real_ else as.double(res$kl),
        status   = if(is.null(res$status)) NA_integer_ else as.integer(res$status),
        error    = if(!is.null(res$error)) res$error else NA_character_
      )
      cat(sprintf("%-20s seed=%d %-12s iters=%4d kl=%.3e status=%s\n",
                  cfg$name, seed, arm$name,
                  results[[length(results)]]$iters,
                  results[[length(results)]]$kl,
                  results[[length(results)]]$status))
    }
  }
}

# Summarize
df_res <- do.call(rbind, lapply(results, as.data.frame))
cat("\n=== SHIP-GATE ANALYSIS ===\n")

unbounded_configs <- c("unbounded_K3", "unbounded_K6")
bounded_configs   <- c("bounded_K6",   "bounded_K9")

# Median iters per arm x regime — explicit loop (names() on list element = NULL)
for (regime_name in c("unbounded", "bounded")) {
  regime_cfgs <- if (regime_name == "unbounded") unbounded_configs else bounded_configs
  sub <- df_res[df_res$config %in% regime_cfgs, ]
  for (arm_name in c("baseline", "fixed_1.2", "fixed_1.4", "adaptive")) {
    med <- median(sub$iters[sub$arm == arm_name], na.rm=TRUE)
    cat(sprintf("  %-10s %-12s median_iters=%.0f\n", regime_name, arm_name, med))
  }
}

divergence_count <- sum(is.na(df_res$iters))
bracket_failures <- sum(df_res$status == -1, na.rm=TRUE)
noconv_flips <- {
  # count (config,seed) pairs where baseline converged but adaptive did not (status=1=NOCONV)
  baseline_ok     <- df_res[df_res$arm=="baseline" & !is.na(df_res$status) & df_res$status==0, c("config","seed")]
  adaptive_noconv <- df_res[df_res$arm=="adaptive" & !is.na(df_res$status) & df_res$status==1, c("config","seed")]
  nrow(merge(baseline_ok, adaptive_noconv))
}

cat(sprintf("\ndivergence_count=%d  bracket_failures=%d  noconv_flips=%d\n",
            divergence_count, bracket_failures, noconv_flips))

# Ship gate numeric summary
ub_sub  <- df_res[df_res$config %in% unbounded_configs, ]
bd_sub  <- df_res[df_res$config %in% bounded_configs, ]
med_ub_baseline  <- median(ub_sub$iters[ub_sub$arm == "baseline"],  na.rm=TRUE)
med_ub_adaptive  <- median(ub_sub$iters[ub_sub$arm == "adaptive"],  na.rm=TRUE)
med_bd_baseline  <- median(bd_sub$iters[bd_sub$arm == "baseline"],  na.rm=TRUE)
med_bd_adaptive  <- median(bd_sub$iters[bd_sub$arm == "adaptive"],  na.rm=TRUE)
med_all_fixed14  <- median(df_res$iters[df_res$arm == "fixed_1.4"], na.rm=TRUE)
med_all_adaptive <- median(df_res$iters[df_res$arm == "adaptive"],  na.rm=TRUE)
pct_ub <- (med_ub_adaptive - med_ub_baseline) / med_ub_baseline * 100
pct_bd <- (med_bd_adaptive - med_bd_baseline) / med_bd_baseline * 100
pct_bf <- (med_all_adaptive - med_all_fixed14) / med_all_fixed14 * 100
cat(sprintf("adaptive_vs_baseline_unbounded_pct = %.1f%%\n", pct_ub))
cat(sprintf("adaptive_vs_baseline_bounded_pct   = %.1f%%\n", pct_bd))
cat(sprintf("adaptive_vs_bestfixed_pct          = %.1f%%\n", pct_bf))

cond1 <- pct_ub < 0
cond2 <- pct_bd <= 5
cond3 <- pct_bf <= 5
cond4 <- divergence_count == 0
cond5 <- bracket_failures == 0
cond6 <- noconv_flips == 0
ship  <- cond1 && cond2 && cond3 && cond4 && cond5 && cond6
cat(sprintf("\nGATE1(ub<0): %s  GATE2(bd<=5%%): %s  GATE3(bf<=5%%): %s\n", cond1, cond2, cond3))
cat(sprintf("GATE4(div=0): %s  GATE5(brf=0): %s  GATE6(noconv=0): %s\n", cond4, cond5, cond6))
cat(sprintf("\nVERDICT: %s\n", if (ship) "SHIP" else "NO-SHIP"))
