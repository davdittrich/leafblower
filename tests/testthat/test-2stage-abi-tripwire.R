test_that("hier_levels_used == 2 when hierarchical enabled (raking/sinkhorn/greenkhorn, refine + exact)", {
  set.seed(99)
  n <- 400L

  # --- refine (mode=0): non-orthogonal DGP is fine ---
  df_refine <- data.frame(
    sex = sample(c("m", "f"), n, replace = TRUE, prob = c(0.5, 0.5)),
    age = sample(c("young", "old"), n, replace = TRUE, prob = c(0.6, 0.4))
  )
  tgt_refine <- list(
    sex = c(m = 0.5, f = 0.5),
    age = c(young = 0.55, old = 0.45)
  )
  hier_refine <- list(
    coarse_mask      = c(1L, 0L),
    min_cell_n       = 5L,
    mode             = 0L,
    outer_tol        = 1e-6,
    outer_iterations = 10L
  )

  # --- exact (mode=1): orthogonal DGP required ---
  # coarse = region (N/S); fine = edu (H/L). Fine levels nest within coarse.
  # In exact mode the fine margin must be orthogonal to the coarse partition.
  # We achieve orthogonality by making edu levels disjoint between regions:
  # region N → only "H"; region S → only "L".
  df_exact <- data.frame(
    region = factor(rep(c("N", "S"), each = n / 2L)),
    edu    = factor(c(rep("H", n / 2L), rep("L", n / 2L)))
  )
  tgt_exact <- list(
    region = c(N = 0.5, S = 0.5),
    edu    = c(H = 0.5, L = 0.5)
  )
  hier_exact <- list(
    coarse_mask      = c(1L, 0L),
    min_cell_n       = 1L,
    mode             = 1L,
    outer_tol        = 1e-6,
    outer_iterations = 5L
  )

  for (method in c("raking", "sinkhorn", "greenkhorn")) {
    r_refine <- suppressWarnings(
      leafblower::harvest(df_refine, tgt_refine, method = method,
                          hierarchical = hier_refine)
    )
    expect_equal(attr(r_refine, "result")$hierarchical_levels_used, 2L,
                 label = paste(method, "refine hier_levels_used"))

    r_exact <- suppressWarnings(
      leafblower::harvest(df_exact, tgt_exact, method = method,
                          hierarchical = hier_exact)
    )
    expect_equal(attr(r_exact, "result")$hierarchical_levels_used, 2L,
                 label = paste(method, "exact hier_levels_used"))
  }
})

test_that("hierarchical = NULL takes early-out: hierarchical_levels_used == 0", {
  set.seed(42)
  n   <- 200L
  df  <- data.frame(
    age = sample(c("young", "old"), n, replace = TRUE, prob = c(0.6, 0.4))
  )
  tgt <- list(age = c(young = 0.5, old = 0.5))

  w <- leafblower::harvest(df, tgt, method = "raking", hierarchical = NULL)
  res <- attr(w, "result")

  # ABI tripwire: hierarchical disabled path must zero-init all diagnostic fields.
  expect_equal(res$hierarchical_levels_used, 0L)
  expect_equal(res$n_cells_total,          0L)
  expect_equal(res$n_cells_skipped,        0L)
  expect_equal(res$n_cells_inherited,      0L)
  expect_equal(res$outer_iterations_used,  0L)
  expect_equal(res$outer_residual_final,   0.0)
})
