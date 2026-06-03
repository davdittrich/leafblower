# e65t.4.5 re-gate: Sinkhorn adaptive omega under ORIS shipping conditions.
# Replicates oris_shipgate_fixture.R: adversarial inverted targets (slow convergence,
# estimator has room) + stepstone_small bounded. NOT the runif near-identity toys
# that starved the estimator in the original (false-negative) 4.5 gate.
#
# Arms: baseline (mode0, omega=1.0), fixed_1.4 (mode0, omega=1.4), adaptive (mode1).
# Gate metric: iteration count (matched config,seed). Verified by orchestrator.
suppressMessages(library(leafblower))

get_iters  <- function(r) { res <- attr(r, "result"); if (is.null(res)) NA_integer_ else as.integer(res$iterations) }
get_status <- function(r) { res <- attr(r, "result"); if (is.null(res)) NA_integer_ else as.integer(res$status) }
get_kl     <- function(r) { res <- attr(r, "result"); if (is.null(res)) NA_real_   else as.double(res$kl) }

arms <- list(
  list(name="baseline",  sk_omega=1.0, mode=0L),
  list(name="fixed_1.4", sk_omega=1.4, mode=0L),
  list(name="adaptive",  sk_omega=1.0, mode=1L)
)

run_arm <- function(df, tgt, arm, mw, minw, maxit) {
  r <- tryCatch(
    harvest(df, tgt, method="sinkhorn", max_iterations=maxit,
            max_weight=mw, min_weight=minw,
            sk_omega=arm$sk_omega, sk_omega_mode_id=arm$mode,
            convergence=list(metric="kl"), accelerate=FALSE),
    error=function(e) NULL)
  if (is.null(r)) return(list(iters=NA_integer_, status=NA_integer_, kl=NA_real_))
  list(iters=get_iters(r), status=get_status(r), kl=get_kl(r))
}

# ---- Adversarial synthetic fixture (ORIS-style inverted targets) ----
make_adv <- function(seed, K=8L, n=5000L) {
  set.seed(seed)
  cols <- lapply(seq_len(K), function(k) {
    p <- if (k <= K/2) c(0.85,0.15) else c(0.15,0.85)
    sample(c("a","b"), n, replace=TRUE, prob=p)
  })
  df <- as.data.frame(cols, stringsAsFactors=FALSE); names(df) <- paste0("m",seq_len(K))
  # invert: target the opposite of the sampled distribution -> slow, estimator engages
  tgt <- lapply(seq_len(K), function(k) {
    if (k <= K/2) c(a=0.15,b=0.85) else c(a=0.85,b=0.15)
  }); names(tgt) <- names(df)
  list(df=df, tgt=tgt)
}

seeds <- c(42L,123L,999L)
cat("=== ADVERSARIAL FIXTURE (n=5000, K=8, inverted targets) ===\n")
adv_rows <- list()
for (regime in list(list(tag="unbounded", mw=1000, minw=0),
                    list(tag="bounded",   mw=2,    minw=0.4))) {
  for (s in seeds) {
    fx <- make_adv(s)
    for (a in arms) {
      o <- run_arm(fx$df, fx$tgt, a, regime$mw, regime$minw, 2000L)
      adv_rows[[length(adv_rows)+1]] <- data.frame(regime=regime$tag, seed=s, arm=a$name,
                                                   iters=o$iters, status=o$status, kl=o$kl)
      cat(sprintf("%-10s seed=%3d %-10s iters=%5s status=%s kl=%.3e\n",
                  regime$tag, s, a$name, o$iters, o$status, o$kl))
    }
  }
}
adv <- do.call(rbind, adv_rows)

# ---- Stepstone bounded fixture (ORIS stepstone_mw5 analog) ----
cat("\n=== STEPSTONE_SMALL (10000 rows, 9 margins, mw=5) ===\n")
ss_rows <- list()
pq <- "tests/testthat/fixtures/stepstone_small.parquet"
tg <- "tests/testthat/fixtures/stepstone_small_targets.rds"
if (file.exists(pq) && file.exists(tg)) {
  suppressMessages(library(arrow))
  ss  <- as.data.frame(arrow::read_parquet(pq))
  sst <- readRDS(tg)
  for (a in arms) {
    o <- run_arm(ss, sst, a, 5, 0, 500L)
    ss_rows[[length(ss_rows)+1]] <- data.frame(arm=a$name, iters=o$iters, status=o$status, kl=o$kl)
    cat(sprintf("stepstone_mw5 %-10s iters=%5s status=%s kl=%.3e\n", a$name, o$iters, o$status, o$kl))
  }
  ss <- do.call(rbind, ss_rows)
} else {
  cat("STEPSTONE FIXTURE MISSING\n"); ss <- NULL
}

# ---- Analysis ----
cat("\n=== GATE ANALYSIS (median iters per regime x arm) ===\n")
med <- function(d, rg, ar) median(d$iters[d$regime==rg & d$arm==ar], na.rm=TRUE)
for (rg in c("unbounded","bounded")) {
  b <- med(adv,rg,"baseline"); f <- med(adv,rg,"fixed_1.4"); a <- med(adv,rg,"adaptive")
  cat(sprintf("  %-10s baseline=%.0f  fixed_1.4=%.0f  adaptive=%.0f  | adaptive_vs_baseline=%+.1f%%  adaptive_vs_fixed=%+.1f%%\n",
              rg, b, f, a, 100*(a-b)/b, 100*(a-f)/f))
}
if (!is.null(ss)) {
  b <- ss$iters[ss$arm=="baseline"]; f <- ss$iters[ss$arm=="fixed_1.4"]; a <- ss$iters[ss$arm=="adaptive"]
  cat(sprintf("  stepstone  baseline=%d  fixed_1.4=%d  adaptive=%d  | adaptive_vs_baseline=%+.1f%%  adaptive_vs_fixed=%+.1f%%\n",
              b, f, a, 100*(a-b)/b, 100*(a-f)/f))
}
div <- sum(is.na(adv$iters)) + (if(!is.null(ss)) sum(is.na(ss$iters)) else 0)
noconv <- sum(adv$status==1, na.rm=TRUE) + (if(!is.null(ss)) sum(ss$status==1, na.rm=TRUE) else 0)
cat(sprintf("\ndivergence/NA=%d  noconv(status=1)=%d\n", div, noconv))
