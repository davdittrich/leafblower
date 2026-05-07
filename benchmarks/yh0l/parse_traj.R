#!/usr/bin/env Rscript
# benchmarks/yh0l/parse_traj.R
# T2: parse traj_R.log + traj_Py.log → trajectory.csv
#
# Available fields per kErrCheckInterval check (every 10 iters + iter 1 + last):
#   iEPPA iter %d: errRp=%.3e          (verbose>=1)  → convergence_metric
#   n_cap_active=%d                    (verbose>=2)   → n_cap_active
#   margin=%d: log10(f) range [%.2f, %.2f]  (verbose>=2, K times)
#
# NOTE: marginal_kl and weight_change are NOT emitted per iter at verbose=2
#   (src/ieppa.cpp verbose=2 only emits errRp + n_cap_active per iter).
#   These columns are present in the CSV but always NA.
# NOTE: max_err is emitted once on the TERMINATE line as max_error=<val>;
#   it is populated only on the final (terminate) iter row.
#
# CSV columns: lang, iter, convergence_metric, max_err, marginal_kl, weight_change, n_cap_active

suppressPackageStartupMessages({
  if (!requireNamespace("readr", quietly = TRUE)) {
    write_csv <- function(x, path) utils::write.csv(x, path, row.names = FALSE)
  } else {
    library(readr)
  }
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

args <- commandArgs(trailingOnly = FALSE)
script_flag <- grep("^--file=", args, value = TRUE)
if (length(script_flag) > 0) {
  SCRIPT_DIR <- dirname(normalizePath(sub("^--file=", "", script_flag[1])))
} else {
  SCRIPT_DIR <- getwd()
}
ROOT <- normalizePath(file.path(SCRIPT_DIR, "..", ".."), mustWork = FALSE)

LOG_R   <- file.path(ROOT, "benchmarks", "yh0l", "traj_R.log")
LOG_PY  <- file.path(ROOT, "benchmarks", "yh0l", "traj_Py.log")
OUT_CSV <- file.path(ROOT, "benchmarks", "yh0l", "trajectory.csv")

# Parse max_error from the TERMINATE line of a log file.
parse_max_err <- function(lines) {
  pat <- "## TERMINATE status=\\S+ \\(\\S+\\) iterations=\\S+ max_error=([0-9eE.+\\-]+|NA)"
  for (ln in rev(lines)) {
    m <- regexec(pat, ln, perl = TRUE)
    parts <- regmatches(ln, m)[[1]]
    if (length(parts) > 1) {
      v <- suppressWarnings(as.numeric(parts[2]))
      return(v)  # NA_real_ if "NA"
    }
  }
  NA_real_
}

parse_log <- function(path, lang) {
  if (!file.exists(path)) stop(sprintf("Log not found: %s", path))
  lines <- readLines(path, warn = FALSE)

  # R log prefix: "[leafblower] iEPPA iter %d: errRp=..."
  # Py log prefix: "iEPPA iter %d: errRp=..."  (no prefix)
  iter_pat <- "(?:^|\\[leafblower\\] )iEPPA iter (\\d+): errRp=([0-9eE.+\\-]+)"
  cap_pat  <- "(?:^|\\[leafblower\\])\\s+n_cap_active=(\\d+)"

  iters   <- integer(0)
  metrics <- numeric(0)
  n_caps  <- integer(0)

  cur_iter   <- NA_integer_
  cur_metric <- NA_real_
  cur_ncap   <- NA_integer_

  flush_cur <- function() {
    if (!is.na(cur_iter)) {
      iters   <<- c(iters,   cur_iter)
      metrics <<- c(metrics, cur_metric)
      n_caps  <<- c(n_caps,  cur_ncap)
    }
  }

  for (ln in lines) {
    m_iter <- regexpr(iter_pat, ln, perl = TRUE)
    if (m_iter > 0) {
      flush_cur()
      parts      <- regmatches(ln, regexec(iter_pat, ln, perl = TRUE))[[1]]
      cur_iter   <- as.integer(parts[2])
      cur_metric <- as.numeric(parts[3])
      cur_ncap   <- NA_integer_
      next
    }
    m_cap <- regexpr(cap_pat, ln, perl = TRUE)
    if (m_cap > 0) {
      parts    <- regmatches(ln, regexec(cap_pat, ln, perl = TRUE))[[1]]
      cur_ncap <- as.integer(parts[2])
    }
  }
  flush_cur()

  if (length(iters) == 0) {
    warning(sprintf("No iter lines found in %s", path))
    return(data.frame(
      lang               = character(),
      iter               = integer(),
      convergence_metric = numeric(),
      max_err            = numeric(),
      marginal_kl        = numeric(),
      weight_change      = numeric(),
      n_cap_active       = integer(),
      stringsAsFactors   = FALSE
    ))
  }

  # max_err: single value from terminate line; populate on last iter row only
  terminate_max_err <- parse_max_err(lines)
  max_err_vec <- rep(NA_real_, length(iters))
  if (!is.na(terminate_max_err)) {
    max_err_vec[length(max_err_vec)] <- terminate_max_err
  }

  data.frame(
    lang               = lang,
    iter               = iters,
    convergence_metric = metrics,
    max_err            = max_err_vec,
    marginal_kl        = NA_real_,
    weight_change      = NA_real_,
    n_cap_active       = n_caps,
    stringsAsFactors   = FALSE
  )
}

cat("[parse_traj] Parsing R log...\n")
df_r  <- parse_log(LOG_R,  "R")
cat(sprintf("[parse_traj] R:  %d iter records\n", nrow(df_r)))

cat("[parse_traj] Parsing Py log...\n")
df_py <- parse_log(LOG_PY, "Py")
cat(sprintf("[parse_traj] Py: %d iter records\n", nrow(df_py)))

traj <- rbind(df_r, df_py)
traj <- traj[order(traj$lang, traj$iter), ]

write_csv(traj, OUT_CSV)
cat(sprintf("[parse_traj] Written: %s (%d rows)\n", OUT_CSV, nrow(traj)))
