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
