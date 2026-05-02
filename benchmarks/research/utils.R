build_design_matrix <- function(df, tgt) {
  n_row <- nrow(df)
  K <- length(tgt)
  
  J_k <- sapply(tgt, length)
  n_col <- sum(J_k)
  
  p <- integer(n_row + 1)
  j <- integer(n_row * K)
  x <- rep(1.0, n_row * K)
  
  p <- seq(0L, n_row * K, by = K)
  
  cum_J_k <- c(0L, cumsum(J_k)[-K])
  
  j_mat <- matrix(0L, nrow = n_row, ncol = K)
  for (k in seq_len(K)) {
    col_idx <- as.integer(df[[names(tgt)[k]]]) - 1L
    j_mat[, k] <- cum_J_k[k] + col_idx
  }
  
  j <- as.vector(t(j_mat))
  
  list(
    p = as.integer(p),
    j = as.integer(j),
    n_row = as.integer(n_row),
    n_col = as.integer(n_col),
    x = as.numeric(x)
  )
}

make_t1_small <- function(seed) {
  set.seed(seed)
  n <- 1000L
  K <- 3L
  ncat <- 3L
  df <- as.data.frame(lapply(seq_len(K), function(k) factor(sample(letters[1:ncat], n, TRUE))))
  names(df) <- paste0("m", seq_len(K))
  d <- rep(1, n)
  tgt <- lapply(df, function(f) prop.table(table(f)))
  A_csr <- build_design_matrix(df, tgt)
  b <- unlist(lapply(tgt, function(t) t * n))
  list(n=n, K=K, df=df, tgt=tgt, A_csr=A_csr, b=as.numeric(b), d=as.numeric(d), 
       lo=rep(0, n), hi=rep(5, n), max_weight=5)
}

make_stepstone_K9 <- function() {
  df <- arrow::read_parquet("benchmarks/stepstone_bench_data.parquet")
  df$uuid <- NULL
  tgt_list <- jsonlite::fromJSON("benchmarks/stepstone_bench_targets.json")
  tgt <- lapply(tgt_list, function(t) { t <- unlist(t); t/sum(t) })
  for(nm in names(tgt)) df[[nm]] <- factor(df[[nm]])
  
  n <- nrow(df)
  K <- length(tgt)
  A_csr <- build_design_matrix(df, tgt)
  b <- unlist(lapply(tgt, function(t) t * n))
  d <- rep(1, n)
  list(n=n, K=K, df=df, tgt=tgt, A_csr=A_csr, b=as.numeric(b), d=as.numeric(d), 
       lo=rep(0, n), hi=rep(5, n), max_weight=5)
}

make_kk1204_K20 <- function(seed) {
  set.seed(seed)
  n <- 1e6L
  K <- 20L
  nj <- 5L
  df <- as.data.frame(lapply(seq_len(K), function(k) factor(sample(letters[seq_len(nj)], n, TRUE))))
  names(df) <- paste0("m", seq_len(K))
  tgt <- lapply(df, function(f) {
    p <- c(0.6, 0.2, 0.1, 0.07, 0.03)
    p <- p / sum(p)
    setNames(p, levels(f))
  })
  
  A_csr <- build_design_matrix(df, tgt)
  b <- unlist(lapply(tgt, function(t) t * n))
  d <- rep(1, n)
  list(n=n, K=K, df=df, tgt=tgt, A_csr=A_csr, b=as.numeric(b), d=as.numeric(d), 
       lo=rep(0, n), hi=rep(3, n), max_weight=3)
}

check_memory <- function(n, K, max_iter) {
  needed <- n * K * 8 + n * 8 * 5 + max_iter * 10 * 8
  avail <- as.numeric(system2("awk", c("'/MemAvailable/ {print $2}'", "/proc/meminfo"), stdout=TRUE)) * 1024
  if (needed > 0.5 * avail) {
    stop(sprintf("insufficient RAM: need %.2f GB; have %.2f GB", needed/1e9, avail/1e9))
  }
}
