test_that("post_distance_within works with intended input",{
  data(fit_data)
  results <- post_distance_within(fit_data[[1]],
                                  pred = FALSE,
                                  n_draws = 100,
                                  comp = "frob")

  expect_type(results, "list")
  expect_true(all(c("comp", "mod_one", "mod_two", "beta", "pcor") %in% colnames(results)))

})

test_that("post_distance_within works with sequential sampling",{
  data(fit_data)
  results <- post_distance_within(fit_data[[1]],
                                  pred = FALSE,
                                  n_draws = 100,
                                  comp = "frob",
                                  sampling_method = "sequential")

  expect_type(results, "list")
  expect_true(all(c("comp", "mod_one", "mod_two", "beta", "pcor") %in% colnames(results)))

})

test_that("post_distance_within works with differenct comps", {
  data(fit_data)
  results <- post_distance_within(fit_data[[1]],
                                  pred = FALSE,
                                  n_draws = 100,
                                  comp = "l1")

  expect_type(results, "list")
  expect_true(all(c("comp", "mod_one", "mod_two", "beta", "pcor") %in% colnames(results)))

  results <- post_distance_within(fit_data[[1]],
                                  pred = FALSE,
                                  n_draws = 100,
                                  comp = "maxdiff")

  expect_type(results, "list")
  expect_true(all(c("comp", "mod_one", "mod_two", "beta", "pcor") %in% colnames(results)))

})

test_that("post_distance_within breaks when input is not a list",{
  expect_error(post_distance_within(matrix(rnorm(100)),
                                    pred = FALSE,
                                    n_draws = 100,
                                    comp = "frob"))
})

test_that("post_distance_within accepts indices",{
  data(fit_data)
  results <- post_distance_within(fit_data[[1]],
                                  pred = FALSE,
                                  n_draws = 100,
                                  comp = "frob",
                                  indices = list(beta = 1:2, pcor = 3:4))

  expect_type(results, "list")
  expect_true(all(c("comp", "mod_one", "mod_two", "beta", "pcor") %in% colnames(results)))

})

test_that("post_distance_within works for posterior predictive lists",{
  data(fit_data)
  # Create fake posterior predictive draws
  pp_list <- list()
  pp_list$fit <- lapply(1:100, function(x) fit_data[[1]])

  results <- post_distance_within(pp_list,
                                  pred = TRUE,
                                  n_draws = 50,
                                  comp = "frob")


  expect_type(results, "list")
  expect_true(all(c("comp", "mod_one", "mod_two", "beta", "pcor") %in% colnames(results)))

})

test_that("post_distance_within works for posterior predictive lists with indices",{
  data(fit_data)
  # Create fake posterior predictive draws
  pp_list <- list()
  pp_list$fit <- lapply(1:100, function(x) fit_data[[1]])

  results <- post_distance_within(pp_list,
                                  pred = TRUE,
                                  n_draws = 50,
                                  comp = "frob",
                                  indices = list(beta = 1:2, pcor = 3:4))

  expect_type(results, "list")
  expect_true(all(c("comp", "mod_one", "mod_two", "beta", "pcor") %in% colnames(results)))

})

test_that("post_distance_within for posterior predictive lists breaks when one model is not a list",{
  data(fit_data)
  # Create fake posterior predictive draws
  pp_list <- list()
  pp_list$fit <- lapply(1:100, function(x) fit_data[[1]])
  pp_list$fit[[1]] <- matrix(rnorm(100))

  expect_error(post_distance_within(pp_list,
                                    pred = TRUE,
                                    n_draws = 50,
                                    comp = "frob",
                                    indices = list(beta = 1:2, pcor = 3:5)))

})
