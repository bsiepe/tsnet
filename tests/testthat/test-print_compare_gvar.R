test_that("print_compare_gvar works for intended input", {
  data(fit_data)
  test_res <-
    compare_gvar(fit_data[[1]], fit_data[[2]], n_draws = 100)
  expect_output(print(test_res))
  expect_no_error(print(test_res))

})


