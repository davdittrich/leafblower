dyn.load("research/leafblower_research.so")
source("benchmarks/research/utils.R")

set.seed(1L)
n <- 200L
K <- 3L
ncat <- 3L
df <- as.data.frame(lapply(seq_len(K), function(k) factor(sample(letters[1:ncat], n, TRUE))))
names(df) <- paste0("m", seq_len(K))
d <- rep(1, n)

tgt <- lapply(df, function(f) prop.table(table(f)))
A_csr <- build_design_matrix(df, tgt)
b <- unlist(lapply(tgt, function(t) t * n))

r_cp <- .Call("cp_solve_R", A_csr, as.numeric(b), as.numeric(d), rep(0, n), rep(5, n), 1000L, FALSE, 1L)

stopifnot(r_cp$status_code == 0L)
max_err <- max(abs(r_cp$weights - d))
stopifnot(max_err < 1e-8)
cat("WU-2 sanity PASS: max_err =", max_err, "\n")
