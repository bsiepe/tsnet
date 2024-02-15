test_that("posterior_plot works with tsnet_fit input",{
  data(ts_data)
  ts_data1 <- ts_data[1:100,1:3]
  stan_fit <- suppressWarnings(stan_gvar(data = ts_data1,
                        n_chains = 1,
                        n_cores = 1))
  expect_no_error(
    post_plot <- posterior_plot(stan_fit,
                                mat = "beta")
  )
  expect_no_error(
    post_plot <- posterior_plot(stan_fit,
                                mat = "pcor")
  )
})

test_that("posterior_plot fails with incorrect CIs",{
  data(ts_data)
  ts_data1 <- ts_data[1:100,1:3]
  stan_fit <- suppressWarnings(stan_gvar(data = ts_data1,
                        n_chains = 1,
                        n_cores = 1))
  expect_error(
    post_plot <- posterior_plot(stan_fit,
                                mat = "beta",
                                ci = c(-0.1, 0.9))
  )

})

test_that("posterior_plot fails with bad column names",{
  data(ts_data)
  ts_data1 <- ts_data[1:100,1:3]
  colnames(ts_data1) <- c("a_1", "b", "c")
  stan_fit <- suppressWarnings(stan_gvar(data = ts_data1,
                        n_chains = 1,
                        n_cores = 1))
  expect_error(
    post_plot <- posterior_plot(stan_fit,
                                mat = "beta")
  )

})

