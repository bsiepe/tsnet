test_that("stan_fit_convert works with intended input",{
  data(ts_data)
  ts_data1 <- ts_data[1:100,1:3]
  stan_fit <- stan_gvar(data = ts_data1,
                   n_chains = 1,
                   n_cores = 1)
  expect_no_error(
    samples <- stan_fit_convert(stan_fit, return_params = c("beta", "pcor", "sigma"))
  )
}
)
test_that("stan_fit returns error with wrong input",{
  mock_fit <- matrix(1:10)
  expect_error(
    samples <- stan_fit_convert(mock_fit, return_params = c("beta", "pcor", "sigma"))
  )
})
