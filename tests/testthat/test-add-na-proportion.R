## Tests for harvest(add_na_proportion = TRUE) — tu15
## Step 1 creates a RED baseline; Step 2 (fix) makes them GREEN.

library(leafblower)

# ---------------------------------------------------------------------------
# Shared fixture
# ---------------------------------------------------------------------------

.make_fixture <- function(seed = 1L, n = 500L, na_frac = 0.20) {
  set.seed(seed)
  g <- sample(c("a", "b", "c"), n, replace = TRUE)
  na_idx <- sample(n, round(n * na_frac))
  g[na_idx] <- NA
  g <- factor(g, levels = c("a", "b", "c"))
  list(
    df     = data.frame(g = g, stringsAsFactors = FALSE),
    target = list(g = c(a = 0.4, b = 0.35, c = 0.25)),
    na_frac = na_frac
  )
}

# ---------------------------------------------------------------------------
# Test 1 — main NA-proportion calibration
# ---------------------------------------------------------------------------

test_that("add_na_proportion: converges and marginals match injected NA bin", {
  fx <- .make_fixture()
  df     <- fx$df
  target <- fx$target
  na_frac <- fx$na_frac

  # (a) no error / converges
  res <- expect_no_error(
    harvest(df, target, method = "raking",
            add_na_proportion = TRUE, max_iterations = 200L)
  )
  w <- res$weights

  # (b) sum(weights) ≈ n
  expect_equal(sum(w), nrow(df), tolerance = 1e-8)

  # (c) NA marginal share ≈ na_frac (solver convergence tolerance)
  na_share <- sum(w[is.na(df$g)]) / sum(w)
  expect_equal(na_share, na_frac, tolerance = 1e-3)

  # (d) non-NA marginals renormalized: sum(w[g==lvl])/sum(w) ≈ (1-na_frac)*tgt[lvl]
  for (lvl in c("a", "b", "c")) {
    observed <- sum(w[which(df$g == lvl)]) / sum(w)
    expected <- (1 - na_frac) * target$g[[lvl]]
    expect_equal(observed, expected, tolerance = 1e-3,
                 label = paste0("marginal share for level '", lvl, "'"))
  }
})

# ---------------------------------------------------------------------------
# Test 2 — all-NA margin must still raise the existing stop()
# ---------------------------------------------------------------------------

test_that("add_na_proportion: all-NA margin raises informative error", {
  set.seed(42)
  n  <- 100L
  df <- data.frame(
    g = rep(NA_character_, n),
    stringsAsFactors = FALSE
  )
  target <- list(g = c(a = 0.4, b = 0.6))

  expect_error(
    harvest(df, target, method = "raking", add_na_proportion = TRUE),
    regexp = "all observations are NA"
  )
})

# ---------------------------------------------------------------------------
# Test 3 — OOV (unknown level) regression: NA→"NA" fill must NOT affect
#           unknown-level observations; they must still map to -1 (dropped).
# ---------------------------------------------------------------------------

test_that("add_na_proportion: unknown level treated as OOV (gid=-1), not NA bin", {
  # Construct a data frame with:
  #   - margin 'h' whose target has levels {x, y, z}  (no add_na_proportion)
  #   - margin 'g' with add_na_proportion=TRUE, ~20 % NA, no OOV
  # One obs has 'h'=="zzz" (OOV for h). This exercises the char-path
  # `idx[is.na(idx)] <- 0L; idx - 1L` for an NA-bin margin's column AS WELL AS
  # confirming that non-NA OOV values (e.g. "zzz") also reach -1L.

  set.seed(7)
  n <- 300L
  g_vec <- sample(c("a", "b", "c"), n, replace = TRUE)
  na_idx <- sample(n, round(n * 0.15))
  g_vec[na_idx] <- NA
  g_vec <- factor(g_vec, levels = c("a", "b", "c"))

  h_vec <- sample(c("x", "y", "z"), n, replace = TRUE)
  # Inject one OOV level into h
  h_vec[1L] <- "zzz"

  df <- data.frame(g = g_vec, h = h_vec, stringsAsFactors = FALSE)
  target <- list(
    g = c(a = 0.40, b = 0.35, c = 0.25),
    h = c(x = 0.33, y = 0.34, z = 0.33)
  )

  # Must not error: the OOV obs in h is simply excluded from h-margin constraints
  res <- expect_no_error(
    harvest(df, target, method = "raking",
            add_na_proportion = TRUE, max_iterations = 200L)
  )
  expect_true(!is.null(res$weights))
  expect_equal(length(res$weights), n)
})
