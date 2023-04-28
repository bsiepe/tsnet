test_that("input check for fit objects", {
  fit_a <- list()
  fit_b <- list()
  expect_error(compare_gvar(fit_a, fit_b))
})
