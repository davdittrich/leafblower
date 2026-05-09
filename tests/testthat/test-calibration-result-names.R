test_that("3qtd: sinkhorn result message contains algorithm name", {
  skip_if_not_installed("leafblower")
  library(leafblower)
  set.seed(42)
  n <- 200
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace=TRUE)),
    b = factor(sample(c("x","y"), n, replace=TRUE))
  )
  target <- list(a=c("1"=1/3,"2"=1/3,"3"=1/3), b=c(x=0.5,y=0.5))
  w <- harvest(data, target, method="sinkhorn", max_iterations=100L, attach_weights=FALSE)
  r <- attr(w, "result")
  expect_true(grepl("sinkhorn", r$message, ignore.case=TRUE),
    info=paste("message was:", r$message))
})

test_that("3qtd: greg result message contains algorithm name", {
  skip_if_not_installed("leafblower")
  library(leafblower)
  set.seed(42)
  n <- 200
  data <- data.frame(a=factor(sample(c("1","2","3"), n, replace=TRUE)))
  target <- list(a=c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- harvest(data, target, method="greg", max_iterations=100L, attach_weights=FALSE)
  r <- attr(w, "result")
  expect_true(grepl("greg", r$message, ignore.case=TRUE),
    info=paste("message was:", r$message))
})

test_that("3qtd: chebyshev result message contains algorithm name", {
  skip_if_not_installed("leafblower")
  library(leafblower)
  set.seed(42)
  n <- 200
  data <- data.frame(a=factor(sample(c("1","2","3"), n, replace=TRUE)))
  target <- list(a=c("1"=1/3,"2"=1/3,"3"=1/3))
  w <- harvest(data, target, method="chebyshev", max_iterations=100L, attach_weights=FALSE)
  r <- attr(w, "result")
  expect_true(grepl("chebyshev", r$message, ignore.case=TRUE),
    info=paste("message was:", r$message))
})
