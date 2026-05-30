test_that("ORIS converges on bound-hitting problem (max_weight=5, skewed sample)", {
  set.seed(42); n <- 10000L
  df <- data.frame(
    age = factor(sample(c("18-34","35-54","55+"), n, replace=TRUE,
                        prob=c(0.60,0.30,0.10))),
    sex = factor(sample(c("M","F"), n, replace=TRUE, prob=c(0.70,0.30))),
    edu = factor(sample(c("HS","Some","BA","Grad"), n, replace=TRUE,
                        prob=c(0.50,0.25,0.15,0.10)))
  )
  tgt <- list(age=c("18-34"=0.33,"35-54"=0.40,"55+"=0.27),
              sex=c(M=0.49,F=0.51),
              edu=c(HS=0.28,Some=0.30,BA=0.27,Grad=0.15))
  result <- harvest(df, tgt, method="oris", max_weight=5, convergence = list(absolute = 1e-6))
  expect_true(max(result$weights) <= 5.0 + 1e-7)
  diag <- diagnose_weights(result, tgt, result$weights)
  expect_true(all(abs(diag$error_weighted) < 1e-6))
})
