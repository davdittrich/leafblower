test_that("A3: SOR auto triggers on oscillatory tight-clamp synthetic", {
  set.seed(31415)
  n <- 5000
  data <- data.frame(
    v1 = factor(sample(c("A","B","C","D"), n, replace = TRUE)),
    v2 = factor(sample(c("X","Y","Z"), n, replace = TRUE)),
    v3 = factor(sample(c("1","2","3","4","5"), n, replace = TRUE)),
    v4 = factor(sample(c("p","q"), n, replace = TRUE)),
    v5 = factor(sample(c("a","b","c","d","e","f"), n, replace = TRUE))
  )
  target <- list(
    v1 = c(A=0.1, B=0.4, C=0.4, D=0.1),
    v2 = c(X=0.5, Y=0.3, Z=0.2),
    v3 = c("1"=0.1,"2"=0.1,"3"=0.4,"4"=0.3,"5"=0.1),
    v4 = c(p=0.7, q=0.3),
    v5 = c(a=0.05,b=0.05,c=0.5,d=0.2,e=0.15,f=0.05)
  )
  w <- leafblower::harvest(data, target, max_weight = 1.5, method = "ieppa",
                           max_iterations = 1000,
                           sor = list(auto = TRUE, omega_min = 0.3),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_lt(result$sor$min_omega, 0.9)
})

test_that("A4: SOR silent on smooth input — no damping", {
  set.seed(202)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("1","2","3"), n, replace = TRUE)),
    b = factor(sample(c("1","2"), n, replace = TRUE))
  )
  target <- list(a = c("1"=1/3,"2"=1/3,"3"=1/3), b = c("1"=0.5,"2"=0.5))
  w <- leafblower::harvest(data, target, max_weight = 10, method = "ieppa",
                           max_iterations = 500,
                           sor = list(auto = TRUE, omega_min = 0.3),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$sor$min_omega, 1.0)
  expect_equal(result$sor$n_damped, 0L)
})
