test_that("all overlays default-off: baseline identical to defaulted", {
  set.seed(1)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("x","y","z"), n, replace = TRUE)),
    b = factor(sample(c("p","q"), n, replace = TRUE))
  )
  target <- list(
    a = c(x = 0.3, y = 0.5, z = 0.2),
    b = c(p = 0.6, q = 0.4)
  )
  baseline <- leafblower::harvest(
    data, target, max_weight = 3, method = "ieppa",
    attach_weights = FALSE
  )
  defaulted <- leafblower::harvest(
    data, target, max_weight = 3, method = "ieppa",
    homotopy_levels = 1,
    scheduler      = "round_robin",
    eta_schedule   = "fixed",
    attach_weights = FALSE
  )
  expect_equal(as.numeric(baseline), as.numeric(defaulted), tolerance = 1e-12)
})

test_that("overlay args run identity-only (scaffolding WU-1)", {
  set.seed(2)
  n <- 2000
  data <- data.frame(
    a = factor(sample(c("x","y","z"), n, replace = TRUE)),
    b = factor(sample(c("p","q"), n, replace = TRUE))
  )
  target <- list(
    a = c(x = 0.3, y = 0.5, z = 0.2),
    b = c(p = 0.6, q = 0.4)
  )
  common <- list(data = data, target = target, max_weight = 3,
                 method = "ieppa", attach_weights = FALSE)
  baseline <- do.call(leafblower::harvest, common)
  a_only <- do.call(leafblower::harvest,
    c(common, list(homotopy_levels = 3, homotopy_start_factor = 5,
                   homotopy_end_factor = 1)))
  b_only <- do.call(leafblower::harvest,
    c(common, list(scheduler = "greedy")))
  e_only <- do.call(leafblower::harvest,
    c(common, list(homotopy_levels = 3, homotopy_start_factor = 5,
                   homotopy_end_factor = 1,
                   eta_schedule = "tang_dynamic",
                   eta_start = 5, eta_end = 1)))
  expect_true(max(as.numeric(a_only)) <= 3 + 1e-10)
  expect_true(max(as.numeric(b_only)) <= 3 + 1e-10)
  expect_true(max(as.numeric(e_only)) <= 3 + 1e-10)
  # Behavioural identity: all overlays no-op in WU-1.
  expect_equal(as.numeric(a_only), as.numeric(baseline), tolerance = 1e-12)
  expect_equal(as.numeric(b_only), as.numeric(baseline), tolerance = 1e-12)
  expect_equal(as.numeric(e_only), as.numeric(baseline), tolerance = 1e-12)
})
