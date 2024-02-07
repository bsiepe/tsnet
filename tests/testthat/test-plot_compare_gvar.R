test_that("plot.compare_gvar works with expected input", {
  data(fit_data)
  test_res <- compare_gvar(fit_data[[1]],
                           fit_data[[2]],
                           n_draws = 100,
                           return_all = TRUE)

  result <- plot.compare_gvar(test_res)

  expect_s3_class(result, "gg")
  expect_no_error(result)

})

test_that("plot.compare_gvar does not work with unexpected input", {
  test_res <- list(
    res_beta = list(mod = "mod_a", null = rnorm(10)),
    res_pcor = list(mod = "mod_b", null = rnorm(10))
  )
  class(test_res) <- "invalid_class"

  # Call the function and check that it throws an error
  expect_error(plot.compare_gvar(test_res))
})

test_that("plot.compare_gvar warns about missing return_all", {
  data(fit_data)
  test_res <- compare_gvar(fit_data[[1]],
                           fit_data[[2]],
                           n_draws = 100,
                           return_all = FALSE)

  expect_error(plot.compare_gvar(test_res))
})

test_that("plot.compare_gvar changes model names", {
  data(fit_data)
  test_res <- compare_gvar(fit_data[[1]],
                                      fit_data[[2]],
                                      n_draws = 100,
                                      return_all = TRUE)

  result <- plot.compare_gvar(test_res,
                              name_a = "Model A",
                              name_b = "Model B")

  expect_s3_class(result, "gg")
  expect_no_error(result)

})

test_that("plot.compare_gvar works with both dec_rules",{
  data(fit_data)
  test_res <- compare_gvar(fit_data[[1]],
                           fit_data[[2]],
                           n_draws = 100,
                           return_all = TRUE,
                           dec_rule = "comb")
  result <- plot.compare_gvar(test_res)
  expect_s3_class(result, "gg")
  expect_no_error(result)
})
