# CR-D16 (leafblower-j7x8.16): chebyshev/greenkhorn/logit now surface the real
# n_bounds_violated count (previously dropped into a local -> stale 0). The three
# result structs gained the field; the has_n_bounds trait then surfaces it.

test_that("chebyshev and greenkhorn surface non-zero n_bounds_violated (CR-D16)", {
  set.seed(3); n <- 400
  df <- data.frame(a  = factor(sample(c("x", "y", "z"), n, TRUE, c(.7, .2, .1))),
                   b  = factor(sample(c("p", "q"), n, TRUE, c(.65, .35))),
                   cc = factor(sample(c("m", "o"), n, TRUE)))
  tg <- list(a = c(x = .34, y = .33, z = .33), b = c(p = .5, q = .5), cc = c(m = .5, o = .5))
  for (m in c("chebyshev", "greenkhorn")) {
    r <- attr(suppressWarnings(
      harvest(df, tg, method = m, min_weight = 0.6, max_weight = 1.5, max_iter = 200)), "result")
    expect_true(is.numeric(r$n_bounds_violated),
                info = sprintf("%s: n_bounds_violated must be surfaced (not NULL)", m))
    expect_gt(r$n_bounds_violated, 0L)   # tight bounds bind on this skewed fixture
  }
})

test_that("logit surfaces n_bounds_violated field (wired; 0 by bounds-construction) (CR-D16)", {
  # logit enforces bounds via the sigmoid link, so per-obs weights stay in [min,max]
  # with uniform base weights -> the count is a genuine 0, but the field must be wired.
  set.seed(4); n <- 300
  df <- data.frame(a = factor(sample(c("x", "y"), n, TRUE)),
                   b = factor(sample(c("p", "q"), n, TRUE)))
  tg <- list(a = c(x = .5, y = .5), b = c(p = .5, q = .5))
  r <- attr(suppressWarnings(
    harvest(df, tg, method = "logit", min_weight = 0.2, max_weight = 3, max_iter = 100)), "result")
  expect_true(is.numeric(r$n_bounds_violated))
  expect_gte(r$n_bounds_violated, 0L)
})
