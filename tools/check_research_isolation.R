# Epic-J WU-1 stub
# Pre-commit + pre-push gate. Exits 1 if package .so contains research symbols.

if (!file.exists("src/leafblower.so")) {
  cat("src/leafblower.so not found; skipping isolation check\n")
  quit(status = 0)
}

syms <- system2("nm", c("-D", "src/leafblower.so"), stdout = TRUE)
forbidden <- c("ipm_solve_R", "ipm_calibrate")
hits <- sapply(forbidden, function(s) any(grepl(s, syms, fixed = TRUE)))

if (any(hits)) {
  stop("Research symbols leaked into package .so: ",
       paste(forbidden[hits], collapse = ", "))
}

cat("research isolation OK\n")
