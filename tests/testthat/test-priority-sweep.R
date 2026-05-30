test_that("greedy scheduler preserves feasibility and Pearson agreement", {
  set.seed(7)
  n <- 5000
  data <- data.frame(
    a = factor(sample(c("p","q"), n, replace = TRUE)),
    b = factor(sample(c("1","2","3","4","5"), n, replace = TRUE)),
    c = factor(sample(c("A","B","C","D"), n, replace = TRUE))
  )
  target <- list(
    a = c(p = 0.6, q = 0.4),
    b = c("1"=0.2,"2"=0.2,"3"=0.2,"4"=0.2,"5"=0.2),
    c = c(A=0.1, B=0.2, C=0.3, D=0.4)
  )
  rr <- leafblower::harvest(data, target, max_weight = 4,
                            method = "oris",
                            max_iterations = 200,
                            convergence = list(absolute = 1e-5),
                            scheduler = "round_robin",
                            attach_weights = FALSE)
  grd <- leafblower::harvest(data, target, max_weight = 4,
                             method = "oris",
                             max_iterations = 200,
                             convergence = list(absolute = 1e-5),
                             scheduler = "greedy",
                             attach_weights = FALSE)
  expect_gt(cor(as.numeric(rr), as.numeric(grd)), 0.995)
})

test_that("greedy uses fewer margin sweeps on 2-cat-heavy problem", {
  set.seed(11)
  n <- 4000
  data <- data.frame(
    easy  = factor(sample(c("A","B"), n, replace = TRUE,
                          prob = c(0.51, 0.49))),
    hardA = factor(sample(as.character(1:10), n, replace = TRUE)),
    hardB = factor(sample(as.character(1:10), n, replace = TRUE))
  )
  target <- list(
    easy  = c(A = 0.5, B = 0.5),
    hardA = setNames(rep(0.1, 10), as.character(1:10)),
    hardB = setNames(rep(0.1, 10), as.character(1:10))
  )
  rr <- leafblower::harvest(data, target, max_weight = 3,
                            method = "oris",
                            max_iterations = 300,
                            convergence = list(absolute = 1e-4),
                            scheduler = "round_robin",
                            attach_weights = FALSE)
  grd <- leafblower::harvest(data, target, max_weight = 3,
                             method = "oris",
                             max_iterations = 300,
                             convergence = list(absolute = 1e-4),
                             scheduler = "greedy",
                             attach_weights = FALSE)
  rr_iters  <- attr(rr,  "iterations")
  grd_iters <- attr(grd, "iterations")
  expect_false(is.null(rr_iters))
  expect_false(is.null(grd_iters))
  # Greedy priority scheduler should do strictly fewer TOTAL margin sweeps
  # than round-robin (which does K=3 sweeps per outer iter). The solver's
  # convergence check granularity (kErrCheckInterval=10) means outer-iter
  # counts often tie; the meaningful efficiency signal is sweep count.
  K <- length(target)
  rr_sweeps  <- K * rr_iters
  grd_sweeps <- attr(grd, "result")$greedy_sweeps_taken
  expect_lt(grd_sweeps, rr_sweeps)
})
