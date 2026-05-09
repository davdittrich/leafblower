# Edge-case validation for design_effect() R wrapper + C++ core.
# Corresponds to validation guards in lbw::design_effect_compute().

test_that("design_effect errors on NA in weights", {
  expect_error(design_effect(c(1, NA, 3)), regexp = "NA|finite", ignore.case = TRUE)
})

test_that("design_effect errors on Inf in weights", {
  expect_error(design_effect(c(1, Inf, 3)), regexp = "NA|finite", ignore.case = TRUE)
})

test_that("design_effect errors when sum(weights) <= 0", {
  expect_error(design_effect(c(0, 0, 0)), regexp = "positive|sum", ignore.case = TRUE)
})

test_that("design_effect 4-arg errors on NA in outcome", {
  w <- c(1, 2, 3); y <- c(10, NA, 30)
  data <- data.frame(g = c("A", "B", "A"), stringsAsFactors = FALSE)
  target <- list(g = c(A = 0.5, B = 0.5))
  expect_error(design_effect(w, outcome = y, data = data, target = target),
               regexp = "NA|finite", ignore.case = TRUE)
})

test_that("design_effect 4-arg errors when outcome length mismatches weights", {
  expect_error(
    design_effect(c(1, 2, 3), outcome = c(10, 20),
                  data = data.frame(g = c("A", "B", "A"), stringsAsFactors = FALSE),
                  target = list(g = c(A = 0.5, B = 0.5))),
    regexp = "length|equal", ignore.case = TRUE
  )
})

test_that("design_effect 4-arg errors when nrow(data) mismatches weights", {
  expect_error(
    design_effect(c(1, 2, 3), outcome = c(10, 20, 30),
                  data = data.frame(g = c("A", "B"), stringsAsFactors = FALSE),
                  target = list(g = c(A = 0.5, B = 0.5))),
    regexp = "nrow|length|equal", ignore.case = TRUE
  )
})

test_that("design_effect 4-arg errors on NA in data column", {
  w <- c(1, 2, 3); y <- c(10, 20, 30)
  data <- data.frame(g = c("A", NA, "B"), stringsAsFactors = FALSE)
  target <- list(g = c(A = 0.5, B = 0.5))
  expect_error(design_effect(w, outcome = y, data = data, target = target),
               regexp = "NA", ignore.case = TRUE)
})

test_that("design_effect 4-arg errors on level not in target", {
  w <- c(1, 2, 3); y <- c(10, 20, 30)
  data <- data.frame(g = c("A", "B", "Z"), stringsAsFactors = FALSE)
  target <- list(g = c(A = 0.5, B = 0.5))
  expect_error(design_effect(w, outcome = y, data = data, target = target),
               regexp = "Z|level|target", ignore.case = TRUE)
})

test_that("design_effect 4-arg collinear margins returns finite deff_H", {
  # Two perfectly-collinear margins. The C++ core uses Gill-Murray-Wright pre-perturbation
  # so dpotrf never fails; result is finite (no crash, no warning from rank_def path).
  set.seed(42)
  n <- 20L; w <- rep(1.0, n); y <- rnorm(n)
  data <- data.frame(
    g1 = rep(c("A", "B"), each = n / 2L),
    g2 = rep(c("A", "B"), each = n / 2L),
    stringsAsFactors = FALSE
  )
  target <- list(g1 = c(A = 0.5, B = 0.5), g2 = c(A = 0.5, B = 0.5))
  d <- design_effect(w, outcome = y, data = data, target = target)
  expect_true(is.finite(d))
  expect_gt(d, 0)
})
