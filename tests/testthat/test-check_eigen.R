test_that("check_eigen works with expected input", {
  # Use the first element of fit_data
  data(fit_data)
  fitobj <- fit_data[[1]]

  # Call the function
  result <- check_eigen(fitobj)

  # Check that the result is a list
  expect_type(result, "list")
})

test_that("check_eigen does not work with unexpected input", {
  # Create a mock fitobj with an invalid class
  fitobj <- list(
    fit = list(
      beta = array(rnorm(1000), dim = c(10, 10, 10)),
      pcors = array(rnorm(1000), dim = c(10, 10, 10))
    )
  )
  class(fitobj) <- "invalid_class"

  # Call the function and check that it throws an error
  expect_error(check_eigen(fitobj))
})

test_that("check_eigen correctly identifies stable models", {
  # Create a mock object with a "stable" matrix
  fitobj <- list(
    beta_mu = matrix(0.05, nrow = 2, ncol = 2)
  )
  class(fitobj) <- "var_estimate"

  # Call the function and check that it identifies the model as stable
  result <- check_eigen(fitobj)
  expect_true(all(Re(result$eigenvalues$eigen_beta_mu)^2 + Im(result$eigenvalues$eigen_beta_mu)^2 < 1))
})

test_that("check_eigen correctly identifies unstable models", {
  # Create a mock fitobj with an unstable model
  fitobj <- list(
    beta_mu = matrix(2, nrow = 10, ncol = 10)
  )
  class(fitobj) <- "var_estimate"

  # Call the function and check that it identifies the model as unstable
  result <- check_eigen(fitobj)
  expect_false(all(Re(result$eigenvalues$eigen_beta_mu)^2 + Im(result$eigenvalues$eigen_beta_mu)^2 < 1))
})

test_that("check_eigen converts stanfit object",{
  data(ts_data)
  example_data <- ts_data[1:100,1:3]
  fit <- stan_gvar(example_data, n_chains = 1)
  expect_no_error(check_eigen(fit))
  expect_output(check_eigen(fit))
})
