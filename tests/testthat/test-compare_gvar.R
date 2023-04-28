test_that("input check for fit objects", {
  fit_a <- list()
  fit_b <- list()
  expect_error(compare_gvar(fit_a, fit_b))
})

test_that("compare_gvar works for var_estimate objects", {
  data("ts_data")
  data1 <- subset(ts_data, id == "ID1")
  data1 <- subset(data1, select = -c(id))
  fit1 <- BGGM::var_estimate(data1)
  data2 <- subset(ts_data, id == "ID2")
  data2 <- subset(data2, select = -c(id))
  fit2 <- BGGM::var_estimate(data2)

  expect_silent(compare_gvar(fit1, fit2))
})
