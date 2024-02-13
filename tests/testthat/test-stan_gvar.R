test_that("stan_gvar works with expected input", {
  # Use time series data
  data(ts_data)
  example_data <- ts_data[1:100,1:3]
  expect_no_error(result <- suppressWarnings(stan_gvar(example_data,
                      n_chains = 1)))

  expect_s4_class(result$stan_fit, "stanfit")

})

test_that("stan_gvar does not work with unexpected input", {
  # Create a mock data frame with non-numeric values
  data <- data.frame(matrix(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j"), nrow = 5))

  expect_error(stan_gvar(data))
})

test_that("stan_gvar works with variational inference", {
  # Use time series data
  data(ts_data)
  example_data <- ts_data[1:100,1:3]
  # Call the function
  expect_no_error(
    result <- suppressWarnings(stan_gvar(example_data,
                      n_chains = 2,
                      method = "variational",
                      iter_sampling = 5000))
    )

  expect_equal(result$stan_fit@stan_args[[1]]$algorithm, "meanfield")
})

test_that("stan_gvar works with both covariance priors", {

  data(ts_data)
  example_data <- ts_data[1:100,1:3]
  expect_no_error(
    result <- suppressWarnings(stan_gvar(example_data,
                      n_chains = 1,
                      cov_prior = "LKJ"))
  )

  expect_equal(result$stan_fit@stanmodel@model_name, "VAR_LKJ")

  expect_no_error(
    result <- suppressWarnings(stan_gvar(example_data,
                      n_chains = 1,
                      cov_prior = "IW"))
    )
  expect_equal(result$stan_fit@stanmodel@model_name, "VAR_wishart")

})

test_that("stan_gvar center_only works",{
  data(ts_data)
  example_data <- ts_data[1:100,1:3]
  expect_no_error(
    result <- suppressWarnings(stan_gvar(example_data,
                      n_chains = 1,
                      center_only = TRUE))
  )

  expect_s4_class(result$stan_fit, "stanfit")
})

test_that("rmv_overnight and beep work for both priors",{
  data(ts_data)
  example_data <- ts_data[1:100,1:3]
  v_beep <- rep(seq(1,4), 25)
  expect_no_error(
    result <- suppressWarnings(stan_gvar(example_data,
                      n_chains = 1,
                      cov_prior = "LKJ",
                      rmv_overnight = TRUE,
                      beep = v_beep))
  )

  expect_s4_class(result$stan_fit, "stanfit")

  expect_no_error(
    result <- suppressWarnings(stan_gvar(example_data,
                      n_chains = 1,
                      cov_prior = "IW",
                      rmv_overnight = TRUE,
                      beep = v_beep))
    )

  expect_s4_class(result$stan_fit, "stanfit")
})


test_that("stan_gvar accepts custom priors",{
  data(ts_data)
  example_data <- ts_data[1:100,1:3]
  K <- ncol(example_data)

  custom_priors <- list(
  prior_Beta_loc = matrix(rep(0.01, K*K), nrow = K, ncol = K),
  prior_Beta_scale = matrix(rep(0.6, K*K), nrow = K, ncol = K),
  prior_Rho_loc = matrix(rep(0.6, K*K), nrow = K, ncol = K),
  prior_Rho_scale = matrix(rep(sqrt(0.6), K*K), nrow = K, ncol = K),
  prior_Rho_marginal = 0.3,
  prior_Eta = 2,
  prior_s = diag(rep(1.01, K))
  )


  expect_no_error(
    result <- suppressWarnings(stan_gvar(example_data,
                      n_chains = 1,
                      priors = custom_priors,
                      cov_prior = "IW"))
  )
  expect_s4_class(result$stan_fit, "stanfit")

  expect_no_error(
    result <- suppressWarnings(stan_gvar(example_data,
                      n_chains = 1,
                      priors = custom_priors,
                      cov_prior = "LKJ"))
  )
  expect_s4_class(result$stan_fit, "stanfit")

})


test_that("stan_gvar accepts additional arguments to rstan", {
  data(ts_data)
  example_data <- ts_data[1:100,1:6]
  result <- suppressWarnings(stan_gvar(example_data,
                      n_chains = 1,
                      control = list(adapt_delta = 0.9,
                                     max_treedepth = 15)
                      ))
  expect_equal(result$stan_fit@stan_args[[1]]$control$adapt_delta, 0.9)
  expect_equal(result$stan_fit@stan_args[[1]]$control$max_treedepth, 15)


})







