test_that("oris+sraa with SOR re-enabled: max_err within 2x of plain", {
  skip_if_not_installed("arrow")
  pq <- testthat::test_path("fixtures/stepstone_small.parquet")
  skip_if_not(file.exists(pq))
  ss <- arrow::read_parquet(pq)
  tgt <- readRDS(testthat::test_path("fixtures/stepstone_small_targets.rds"))
  for (nm in names(tgt)) ss[[nm]] <- factor(ss[[nm]])
  r_plain <- suppressWarnings(harvest(ss, tgt, method = "oris",
                                      max_iterations = 200L,
                                      accelerate = FALSE))
  r_sraa  <- suppressWarnings(harvest(ss, tgt, method = "oris",
                                      max_iterations = 200L,
                                      accelerate = TRUE))
  me_plain <- attr(r_plain, "result")$max_error
  me_sraa  <- attr(r_sraa,  "result")$max_error
  expect_lte(me_sraa, 2 * me_plain + 1e-4,
             label = sprintf("sraa max_err %.4g <= 2x plain %.4g",
                             me_sraa, me_plain))
  expect_equal(sum(r_sraa$weights), nrow(ss), tolerance = 1e-6)
})
