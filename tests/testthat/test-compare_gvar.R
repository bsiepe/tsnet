test_that("input check for fit objects", {
  fit_a <- list()
  fit_b <- list()
  expect_error(compare_gvar(fit_a, fit_b))
})

test_that("compare_gvar fails if one model is incorrect", {
  data(fit_data)
  expect_error(compare_gvar(fit_data[[1]], matrix(rnorm(10)), n_draws = 100))
})


test_that("compare_gvar works for two fits", {
  data(fit_data)
  test_res <-
    compare_gvar(fit_data[[1]], fit_data[[2]], n_draws = 100)
  expect_s3_class(test_res, "compare_gvar")
  expect_no_error(test_res)

})



test_that("compare_gvar returns zero difference for identical models", {
  data(fit_data)
  test_res <-
    compare_gvar(fit_data[[1]], fit_data[[1]], n_draws = 100)
  expect_equal(test_res$emp_beta, 0)
  expect_equal(test_res$emp_pcor, 0)
})


test_that("compare_gvar returns non-zero difference for different models",
          {
            data(fit_data)
            test_res <-
              compare_gvar(fit_data[[1]], fit_data[[2]], n_draws = 100)
            expect_gt(test_res$emp_beta, 0)
            expect_gt(test_res$emp_pcor, 0)
          })

test_that("compare_gvar fails for incorrect cutoff", {
  data(fit_data)
  expect_error(compare_gvar(fit_data[[1]],
                            fit_data[[2]],
                            cutoff = 101,
                            n_draws = 100))
  expect_error(compare_gvar(
    fit_data[[1]],
    fit_data[[2]],
    cutoff = -0.1,
    n_draws = 100
  ))
})

test_that("compare_gvar warns for low cutoff", {
  data(fit_data)
  expect_warning(compare_gvar(fit_data[[1]],
                              fit_data[[2]],
                              cutoff = 0.1,
                              n_draws = 100))
})

test_that("compare_gvar fails for incorrect indices", {
  data(fit_data)
  expect_error(compare_gvar(fit_data[[1]],
                            fit_data[[2]],
                            n_draws = 100,
                            indices = 1:10))
  expect_error(compare_gvar(fit_data[[1]],
                            fit_data[[2]],
                            n_draws = 100,
                            indices = list(theta = 1:10)))
  expect_error(compare_gvar(fit_data[[1]],
                            fit_data[[2]],
                            n_draws = 100,
                            indices = list(theta = 1:10, beta = 1:10)))
  expect_error(compare_gvar(fit_data[[1]],
                            fit_data[[2]],
                            n_draws = 100,
                            indices = list(beta = "wrong", pcor = 1:10)))
  expect_error(compare_gvar(fit_data[[1]],
                            fit_data[[2]],
                            n_draws = 100,
                            indices = list(beta = 1:10, pcor = "wrong")))
})

test_that("compare_gvar breaks for wrong test arguments", {
  data(fit_data)
  expect_error(compare_gvar(fit_data[[1]],
                            fit_data[[2]],
                            n_draws = 100,
                            comp = "wrong"))
  expect_error(compare_gvar(fit_data[[1]],
                            fit_data[[2]],
                            n_draws = 100,
                            dec_rule = "wrong"))
})


test_that("compare_gvar works for correct indices", {
  data(fit_data)
  test_res <-
    compare_gvar(fit_data[[1]],
                 fit_data[[2]],
                 n_draws = 100,
                 indices = list(beta = 1:2, pcor = 1:2))
  expect_s3_class(test_res, "compare_gvar")
  expect_no_error(test_res)
})

test_that("compare_gvar works with different comps", {
  data(fit_data)
  test_res <-
    compare_gvar(fit_data[[1]],
                 fit_data[[2]],
                 n_draws = 100,
                 comp = "l1")
  expect_s3_class(test_res, "compare_gvar")
  expect_no_error(test_res)

  test_res <-
    compare_gvar(fit_data[[1]],
                 fit_data[[2]],
                 n_draws = 100,
                 comp = "maxdiff")
  expect_s3_class(test_res, "compare_gvar")
  expect_no_error(test_res)
})

test_that("compare_gvar works with different dec_rule", {
  data(fit_data)
  test_res <-
    compare_gvar(fit_data[[1]],
                 fit_data[[2]],
                 n_draws = 100,
                 dec_rule = "comb")
  expect_s3_class(test_res, "compare_gvar")
  expect_no_error(test_res)

})


test_that("compare_gvar return_all works",{
  data(fit_data)
  test_res <-
    compare_gvar(fit_data[[1]],
                 fit_data[[2]],
                 n_draws = 100,
                 return_all = TRUE)
  expect_s3_class(test_res, "compare_gvar")
  expect_no_error(test_res)
  expect_true(all(c("res_beta", "res_pcor") %in% names(test_res)))
})

test_that("compare_gvar converts stanfit",{
  data(ts_data)
  example_data <- ts_data[1:100,1:3]
  fit <- stan_gvar(example_data, n_chains = 1)
  test_res <-
    compare_gvar(fit,
                 fit,
                 n_draws = 100)
  expect_s3_class(test_res, "compare_gvar")
  expect_no_error(test_res)
})

test_that("compare_gvar returns expected distances",{
  data(fit_data)
  test_res <-
    compare_gvar(fit_data[[1]],
                 fit_data[[2]],
                 n_draws = 100,
                 comp = "frob")

  # calculate distances between posterior means
  beta_diff <- norm(fit_data[[1]]$beta_mu - fit_data[[2]]$beta_mu,
                    type = "F")
  ut <- function(x) {
      matrix(x[upper.tri(x, diag = FALSE)])
  }
  pcor_diff <- norm(ut(fit_data[[1]]$pcor_mu) - ut(fit_data[[2]]$pcor_mu),
                    type = "F")

  expect_equal(test_res$emp_beta, beta_diff)
  expect_equal(test_res$emp_pcor, pcor_diff)


})
