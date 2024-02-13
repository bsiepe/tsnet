test_that("print.compare_gvar works for intended input", {
  data(fit_data)
  test_res <-
    compare_gvar(fit_data[[1]], fit_data[[2]], n_draws = 100)
  expect_output(print.compare_gvar(test_res))
  expect_no_error(print.compare_gvar(test_res))

})

test_that("print.copmare_gvar gives error for unintended input", {
  data(fit_data)
  expect_error(print.compare_gvar(fit_data[[1]]))
})


test_that("print.tsnet_fit works for intended input", {
  data(ts_data)
  example_data <- ts_data[1:100,1:3]
  expect_no_error(result_mcmc <- suppressWarnings(stan_gvar(example_data,
                                                       n_chains = 1,
                                                       method = "sampling")))
  expect_output(print(result_mcmc))
  expect_no_error(print(result_mcmc))

  expect_no_error(result_vb <- suppressWarnings(stan_gvar(example_data,
                                                            n_chains = 1,
                                                            method = "variational")))
  expect_output(print(result_vb))
  expect_no_error(print(result_vb))


})
