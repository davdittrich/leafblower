test_that("A3: SOR auto triggers on oscillatory tight-clamp synthetic", {
  # Uses mode_id=2 (spectral/free-subspace θ₂) explicitly — that mode damps on
  # oscillatory inputs. Default is mode_id=1 (fixed) after e18t.5 NO-SHIP decision.
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
  w <- leafblower::harvest(data, target, max_weight = 1.5, method = "oris",
                           max_iterations = 1000,
                           convergence = list(absolute = 1e-6),
                           sor = list(auto = TRUE, omega_min = 0.3,
                                      omega_mode_id = 2L),
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
  w <- leafblower::harvest(data, target, max_weight = 10, method = "oris",
                           max_iterations = 500,
                           sor = list(auto = TRUE, omega_min = 0.3),
                           attach_weights = FALSE)
  result <- attr(w, "result")
  expect_equal(result$sor$min_omega, 1.0)
  expect_equal(result$sor$n_damped, 0L)
})

# e18t.5 ship-gate: NO-SHIP (2026-06-02)
# mode_id=2 (spectral/free-subspace θ₂) over-relaxes the stepstone fixture to
# a fixed-point stall (status=4, 500 iters vs 140 for fixed). Condition 2 FAIL.
# mode_id=1 (fixed omega_max) remains the default per the NO-SHIP decision.
test_that("E18T.5: mode_id=1 default converges T2 slow-unconstrained in <= 500 iters", {
  skip_on_cran()
  set.seed(20260531L)
  n <- 5000L
  df <- data.frame(
    m1 = sample(c("a","b"), n, replace = TRUE, prob = c(0.85, 0.15)),
    m2 = sample(c("a","b"), n, replace = TRUE, prob = c(0.85, 0.15)),
    m3 = sample(c("a","b"), n, replace = TRUE, prob = c(0.85, 0.15)),
    m4 = sample(c("a","b"), n, replace = TRUE, prob = c(0.85, 0.15)),
    m5 = sample(c("a","b"), n, replace = TRUE, prob = c(0.15, 0.85)),
    m6 = sample(c("a","b"), n, replace = TRUE, prob = c(0.15, 0.85)),
    m7 = sample(c("a","b"), n, replace = TRUE, prob = c(0.15, 0.85)),
    m8 = sample(c("a","b"), n, replace = TRUE, prob = c(0.15, 0.85))
  )
  tgt <- list(
    m1 = c(a = 0.15, b = 0.85), m2 = c(a = 0.15, b = 0.85),
    m3 = c(a = 0.15, b = 0.85), m4 = c(a = 0.15, b = 0.85),
    m5 = c(a = 0.85, b = 0.15), m6 = c(a = 0.85, b = 0.15),
    m7 = c(a = 0.85, b = 0.15), m8 = c(a = 0.85, b = 0.15)
  )
  # mode_id=1 default: fixed omega_max
  w_fix <- suppressWarnings(
    leafblower::harvest(df, tgt, method = "oris",
                        sor = list(auto = TRUE, omega_mode_id = 1L),
                        max_weight = 1000, min_weight = 0, max_iterations = 2000)
  )
  r_fix <- attr(w_fix, "result")
  # pre-registered baseline: 350 iters
  expect_lte(r_fix$iterations, 500L)
  expect_true(r_fix$status %in% c(0L, 5L))

  # mode_id=2 (spectral) is faster on THIS fixture but NO-SHIP on stepstone
  w_spec <- suppressWarnings(
    leafblower::harvest(df, tgt, method = "oris",
                        sor = list(auto = TRUE, omega_mode_id = 2L),
                        max_weight = 1000, min_weight = 0, max_iterations = 2000)
  )
  r_spec <- attr(w_spec, "result")
  # spectral is faster on T2 (280 < 350)
  expect_lt(r_spec$iterations, r_fix$iterations)
})
