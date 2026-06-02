#!/usr/bin/env Rscript
# oris_baseline_snapshot.R
# Convergence-safety baseline: run canonical ORIS problems under fixed-omega
# (omega_mode_id=1L) and no-SOR, record status+iters+max_error.
# Output: tests/testthat/fixtures/oris_fixed_omega_baseline.json
# Ticket: leafblower-e18t.2

Sys.setenv(OMP_NUM_THREADS    = "1",
           OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS    = "1")

library(leafblower)
library(jsonlite)

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------

run_case <- function(df, tgt, extra = list()) {
  sor_fixed <- list(auto = TRUE, omega_mode_id = 1L)
  base_args <- c(
    list(data = df, target = tgt, method = "oris", max_iterations = 500L),
    extra
  )

  run_one <- function(sor_arg) {
    args <- c(base_args, list(sor = sor_arg))
    res <- tryCatch(
      suppressWarnings(do.call(harvest, args)),
      error = function(e) e
    )
    if (inherits(res, "error")) {
      # status=2 (INFEAS) manifests as an R error; treat as status 2.
      list(status = 2L, iters = NA_integer_, max_error = NA_real_)
    } else {
      info <- attr(res, "result")
      list(
        status    = as.integer(info$status),
        iters     = as.integer(info$iterations),
        max_error = as.double(info$max_error)
      )
    }
  }

  list(
    fixed_omega = run_one(sor_fixed),
    no_sor      = run_one(NULL)
  )
}

converged <- function(status) status %in% c(0L, 5L)

# --------------------------------------------------------------------------
# Fixture 1: simple 2-margin balanced (should converge fast)
# --------------------------------------------------------------------------
set.seed(1L)
n1 <- 2000L
df1 <- data.frame(
  a = sample(c("A","B"), n1, replace = TRUE),
  b = sample(c("X","Y"), n1, replace = TRUE)
)
tgt1 <- list(a = c(A = 0.5, B = 0.5),
             b = c(X = 0.5, Y = 0.5))
fx1 <- run_case(df1, tgt1)
cat("fx1 (simple_2margin_balanced) done\n")

# --------------------------------------------------------------------------
# Fixture 2: simple 2-margin skewed (medium difficulty)
# --------------------------------------------------------------------------
set.seed(2L)
n2 <- 3000L
df2 <- data.frame(
  a = sample(c("A","B"), n2, replace = TRUE, prob = c(0.8, 0.2)),
  b = sample(c("X","Y"), n2, replace = TRUE, prob = c(0.7, 0.3))
)
tgt2 <- list(a = c(A = 0.2, B = 0.8),
             b = c(X = 0.3, Y = 0.7))
fx2 <- run_case(df2, tgt2)
cat("fx2 (simple_2margin_skewed) done\n")

# --------------------------------------------------------------------------
# Fixture 3: T2 ship-gate (n=5000, K=8, seed=20260531, max_weight=1000)
# --------------------------------------------------------------------------
set.seed(20260531L)
n3 <- 5000L
df3 <- data.frame(
  m1 = sample(c("a","b"), n3, replace = TRUE, prob = c(0.85, 0.15)),
  m2 = sample(c("a","b"), n3, replace = TRUE, prob = c(0.85, 0.15)),
  m3 = sample(c("a","b"), n3, replace = TRUE, prob = c(0.85, 0.15)),
  m4 = sample(c("a","b"), n3, replace = TRUE, prob = c(0.85, 0.15)),
  m5 = sample(c("a","b"), n3, replace = TRUE, prob = c(0.15, 0.85)),
  m6 = sample(c("a","b"), n3, replace = TRUE, prob = c(0.15, 0.85)),
  m7 = sample(c("a","b"), n3, replace = TRUE, prob = c(0.15, 0.85)),
  m8 = sample(c("a","b"), n3, replace = TRUE, prob = c(0.15, 0.85))
)
tgt3 <- list(
  m1 = c(a = 0.15, b = 0.85), m2 = c(a = 0.15, b = 0.85),
  m3 = c(a = 0.15, b = 0.85), m4 = c(a = 0.15, b = 0.85),
  m5 = c(a = 0.85, b = 0.15), m6 = c(a = 0.85, b = 0.15),
  m7 = c(a = 0.85, b = 0.15), m8 = c(a = 0.85, b = 0.15)
)
fx3 <- run_case(df3, tgt3, list(max_weight = 1000))
cat("fx3 (t2_shipgate_k8_n5000) done\n")

# --------------------------------------------------------------------------
# Fixture 4: stepstone_small (n=10000, K=9 parquet fixture)
# --------------------------------------------------------------------------
pq_path  <- file.path(dirname(dirname(sys.frame(0)$ofile %||%
                                "tests/testthat/fixtures/x")),
                      "testthat/fixtures/stepstone_small.parquet")
# Resolve relative to this script's location
script_dir <- tryCatch(
  normalizePath(dirname(sys.frame(sys.nframe())$ofile)),
  error = function(e) normalizePath("tests/testthat/fixtures")
)
pq_path  <- file.path(script_dir, "stepstone_small.parquet")
rds_path <- file.path(script_dir, "stepstone_small_targets.rds")

fx4 <- if (file.exists(pq_path) && requireNamespace("arrow", quietly = TRUE)) {
  ss  <- as.data.frame(arrow::read_parquet(pq_path))
  tgt <- readRDS(rds_path)
  for (nm in names(tgt)) ss[[nm]] <- factor(ss[[nm]])
  run_case(ss, tgt)
} else {
  cat("fx4 SKIPPED (arrow or parquet not available)\n")
  NULL
}
if (!is.null(fx4)) cat("fx4 (stepstone_small) done\n")

# --------------------------------------------------------------------------
# Fixture 5: tight-bounds problem (max_weight=2, 3 margins)
# --------------------------------------------------------------------------
set.seed(31415L)
n5 <- 5000L
df5 <- data.frame(
  v1 = factor(sample(c("A","B","C","D"), n5, replace = TRUE)),
  v2 = factor(sample(c("X","Y","Z"),     n5, replace = TRUE)),
  v3 = factor(sample(c("p","q"),         n5, replace = TRUE))
)
tgt5 <- list(
  v1 = c(A = 0.1, B = 0.4, C = 0.4, D = 0.1),
  v2 = c(X = 0.5, Y = 0.3, Z = 0.2),
  v3 = c(p = 0.7, q = 0.3)
)
fx5 <- run_case(df5, tgt5, list(max_weight = 2.0))
cat("fx5 (tight_bounds_mw2) done\n")

# --------------------------------------------------------------------------
# Fixture 6: persistent-infeas (empty target cell → status 2)
# --------------------------------------------------------------------------
set.seed(99L)
n6 <- 500L
df6 <- data.frame(
  a = sample(letters[1:2], n6, replace = TRUE),
  b = sample(letters[1:2], n6, replace = TRUE)
)
# 3rd category 'c' has zero observations → structurally infeasible
tgt6 <- list(
  a = c(a = 0.4, b = 0.3, c = 0.3),
  b = c(a = 0.5, b = 0.5)
)
fx6 <- run_case(df6, tgt6, list(max_weight = 5, min_weight = 0))
cat("fx6 (persistent_infeas) done\n")

# --------------------------------------------------------------------------
# Fixture 7: multi-margin oscillatory (SOR test shape, K=5)
# --------------------------------------------------------------------------
set.seed(42L)
n7 <- 2000L
df7 <- data.frame(
  a = sample(letters[1:3], n7, replace = TRUE, prob = c(0.5, 0.3, 0.2)),
  b = sample(letters[1:3], n7, replace = TRUE, prob = c(0.2, 0.5, 0.3)),
  c = sample(letters[1:4], n7, replace = TRUE, prob = c(0.3, 0.3, 0.2, 0.2)),
  d = sample(letters[1:3], n7, replace = TRUE, prob = c(0.4, 0.3, 0.3)),
  e = sample(letters[1:3], n7, replace = TRUE, prob = c(0.25, 0.25, 0.5))
)
tgt7 <- list(
  a = c(a = 0.40, b = 0.35, c = 0.25),
  b = c(a = 0.30, b = 0.40, c = 0.30),
  c = c(a = 0.25, b = 0.25, c = 0.25, d = 0.25),
  d = c(a = 0.33, b = 0.33, c = 0.34),
  e = c(a = 0.30, b = 0.35, c = 0.35)
)
fx7 <- run_case(df7, tgt7, list(max_weight = 5, min_weight = 0))
cat("fx7 (k5_overlapping_margins) done\n")

# --------------------------------------------------------------------------
# Assemble output
# --------------------------------------------------------------------------

build_entry <- function(name, fx) {
  if (is.null(fx)) return(NULL)
  list(
    name        = name,
    fixed_omega = fx$fixed_omega,
    no_sor      = fx$no_sor
  )
}

all_fixtures <- list(
  build_entry("simple_2margin_balanced",  fx1),
  build_entry("simple_2margin_skewed",    fx2),
  build_entry("t2_shipgate_k8_n5000",    fx3),
  build_entry("stepstone_small",          fx4),
  build_entry("tight_bounds_mw2",         fx5),
  build_entry("persistent_infeas",        fx6),
  build_entry("k5_overlapping_margins",   fx7)
)
all_fixtures <- Filter(Negate(is.null), all_fixtures)

n_conv_fixed <- sum(sapply(all_fixtures,
                            function(e) converged(e$fixed_omega$status)))
n_conv_nosor  <- sum(sapply(all_fixtures,
                            function(e) converged(e$no_sor$status)))

out <- list(
  generated_at             = format(Sys.Date(), "%Y-%m-%d"),
  n_fixtures               = length(all_fixtures),
  n_converged_fixed_omega  = n_conv_fixed,
  n_converged_no_sor       = n_conv_nosor,
  fixtures                 = all_fixtures
)

out_path <- file.path(script_dir, "oris_fixed_omega_baseline.json")
writeLines(jsonlite::toJSON(out, pretty = TRUE, auto_unbox = TRUE, na = "null"),
           con = out_path)

cat(sprintf("\nBaseline written to: %s\n", out_path))
cat(sprintf("Fixtures: %d | converged fixed-omega: %d | converged no-SOR: %d\n",
            length(all_fixtures), n_conv_fixed, n_conv_nosor))

# Print table
cat("\nPer-fixture summary:\n")
fmt_iters <- function(x) if (length(x) == 0L || is.na(x)) "NA" else as.character(x)
fmt_err   <- function(x) if (length(x) == 0L || is.na(x)) "NA" else formatC(x, format="e", digits=3)
for (e in all_fixtures) {
  cat(sprintf("  %-32s  fixed_omega: status=%d iters=%-4s max_err=%s  |  no_sor: status=%d iters=%-4s max_err=%s\n",
              e$name,
              e$fixed_omega$status,
              fmt_iters(e$fixed_omega$iters),
              fmt_err(e$fixed_omega$max_error),
              e$no_sor$status,
              fmt_iters(e$no_sor$iters),
              fmt_err(e$no_sor$max_error)
  ))
}
