test_that("get_centrality returns expected output", {
  # Use the first element of fit_data
  data(fit_data)
  fitobj <- fit_data[[1]]

  # Call the function
  result <- get_centrality(fitobj)

  # Check that the result is a list
  expect_type(result, "list")

  # Check that the list has the expected names
  expect_equal(names(result), c("instrength", "outstrength", "strength", "density_beta", "density_pcor"))

})

test_that("get_centrality throws an error with invalid input", {
  # Create a mock fitobj with an invalid class
  fitobj <- list(
    Y = matrix(rnorm(100), ncol = 10),
    fit = list(
      beta = array(rnorm(1000), dim = c(10, 10, 10)),
      pcors = array(rnorm(1000), dim = c(10, 10, 10))
    )
  )
  class(fitobj) <- "invalid_class"

  # Call the function and check that it throws an error
  expect_error(get_centrality(fitobj))
})

test_that("get_centrality converts stanfit",{
  data(ts_data)
  example_data <- ts_data[1:100,1:3]
  fit <- stan_gvar(example_data, n_chains = 1)
  cent <- get_centrality(fit)
  expect_type(cent, "list")
  expect_no_error(cent)
})


test_that("plot_centrality returns expected output", {
  # Use the first element of fit_data
  data(fit_data)
  fitobj <- fit_data[[1]]

  # Call the function
  cent <- get_centrality(fitobj)
  result <- plot_centrality(cent)

  # check that the result is a ggplot object
  expect_true(is_ggplot(result))

})




test_that("plot_centrality throws an error with invalid input", {
  # Create a mock fitobj with an invalid class
  fitobj <- list(
    Y = matrix(rnorm(100), ncol = 10),
    fit = list(
      beta = array(rnorm(1000), dim = c(10, 10, 10)),
      pcors = array(rnorm(1000), dim = c(10, 10, 10))
    )
  )
  class(fitobj) <- "invalid_class"

  # Call the function and check that it throws an error
  expect_error(plot_centrality(fitobj))
})

test_that("plot_centrality throws an error with wrong CIs", {
  # Use the first element of fit_data
  data(fit_data)
  fitobj <- fit_data[[1]]

  # Call the function
  cent <- get_centrality(fitobj)
  result <- plot_centrality(cent)

  # Call the function and check that it throws an error
  expect_error(plot_centrality(cent, ci = 1,1))}
)

test_that("plot_centrality density plot works", {
  # Use the first element of fit_data
  data(fit_data)
  fitobj <- fit_data[[1]]

  # Call the function
  cent <- get_centrality(fitobj)
  result <- plot_centrality(cent, plot_type = "density")

  # Check that the result is a ggplot object
  expect_no_error(print(result))

})


