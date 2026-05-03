library(leafblower)

test_that("WH-e: newton_tsvd_ratio default 1e-8 matches pre-existing behavior", {
  set.seed(11L); n <- 1000L
  df <- data.frame(
    a = factor(sample(c("x","y","z"), n, TRUE)),
    b = factor(sample(c("M","F"), n, TRUE)),
    c = factor(sample(c("Y","O","S"), n, TRUE))
  )
  tgt <- list(
    a = c(x=0.4, y=0.35, z=0.25),
    b = c(M=0.50, F=0.50),
    c = c(Y=0.45, O=0.40, S=0.15)
  )
  r_default <- harvest(df, tgt, method="newton_kl", attach_weights=FALSE)
  r_explicit <- harvest(df, tgt, method="newton_kl",
                        newton_tsvd_ratio=1e-8, attach_weights=FALSE)
  expect_equal(attr(r_default, "result")$max_error,
               attr(r_explicit, "result")$max_error,
               tolerance=1e-15,
               label="WH-e: default 1e-8 == explicit 1e-8")
})

test_that("WH-e: newton_tsvd_ratio higher value drops more directions", {
  set.seed(11L); n <- 1000L
  df <- data.frame(
    a = factor(sample(c("x","y","z"), n, TRUE)),
    b = factor(sample(c("M","F"), n, TRUE)),
    c = factor(sample(c("Y","O","S"), n, TRUE))
  )
  tgt <- list(
    a = c(x=0.4, y=0.35, z=0.25),
    b = c(M=0.50, F=0.50),
    c = c(Y=0.45, O=0.40, S=0.15)
  )
  r_tight <- harvest(df, tgt, method="newton_kl",
                     newton_tsvd_ratio=1e-12, attach_weights=FALSE)
  r_loose <- harvest(df, tgt, method="newton_kl",
                     newton_tsvd_ratio=1e-2, attach_weights=FALSE)
  n_proj_tight <- attr(r_tight, "result")$n_projected_dims
  n_proj_loose <- attr(r_loose, "result")$n_projected_dims
  expect_gte(n_proj_loose, n_proj_tight,
             label=sprintf("WH-e: loose ratio (%d) >= tight ratio (%d) projected dims",
                           n_proj_loose, n_proj_tight))
})
