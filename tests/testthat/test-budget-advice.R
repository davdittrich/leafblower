library(leafblower)

# CR-F3 (dtkn.3): the stall/BUDGET advice string must name only method values
# that map_method actually accepts. It previously suggested method='oris+accel',
# which map_method rejects via match_exact -> actively misleading users.
test_that("stall advice names only valid map_method values (dtkn.3)", {
  src <- paste(deparse(body(harvest)), collapse = "\n")
  expect_false(grepl("oris\\+accel", src))          # rejected token gone
  expect_true(grepl("method='oris'", src))           # valid replacement present
  expect_identical(leafblower:::map_method("oris"), "oris")
  expect_identical(leafblower:::map_method("newton_kl"), "newton_kl")
  expect_error(leafblower:::map_method("oris+accel"), "must be exactly one of")
})
