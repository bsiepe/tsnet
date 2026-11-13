test_that("gvar_priors returns the defaults of stan_gvar", {
  K <- 4
  # defaults as hard-coded in stan_gvar() up to version 0.2.0
  old_defaults <- list(
    prior_Beta_loc = matrix(rep(0, K * K), nrow = K, ncol = K),
    prior_Beta_scale = matrix(rep(0.5, K * K), nrow = K, ncol = K),
    prior_Rho_loc = matrix(rep(0.5, K * K), nrow = K, ncol = K),
    prior_Rho_scale = matrix(rep(sqrt(0.5), K * K), nrow = K, ncol = K),
    prior_Rho_marginal = 0.25,
    prior_S = diag(K),
    prior_Eta = 1
  )

  expect_identical(gvar_priors(K), old_defaults)
})

test_that("gvar_priors applies overrides", {
  loc <- matrix(c(0.3, 0.2, 0, 0.1), nrow = 2, byrow = TRUE)
  priors <- gvar_priors(2, prior_Beta_loc = loc, prior_Eta = 2)

  expect_identical(priors$prior_Beta_loc, loc)
  expect_identical(priors$prior_Eta, 2)
  expect_identical(priors$prior_Beta_scale, gvar_priors(2)$prior_Beta_scale)
})

test_that("gvar_priors rejects invalid input", {
  expect_error(gvar_priors(3, prior_s = diag(3)), "Unknown prior")
  expect_error(gvar_priors(3, diag(3)), "named")
  expect_error(gvar_priors(3, prior_Beta_loc = diag(2)), "3 x 3")
  expect_error(gvar_priors(3, prior_Beta_scale = diag(3)), "positive")
  expect_error(gvar_priors(3, prior_Rho_loc = diag(3)), "between 0 and 1")
  expect_error(gvar_priors(3, prior_Rho_marginal = 0.6), "at most 0.5")
  expect_error(gvar_priors(3, prior_Eta = 1.5), "integer")
  expect_error(gvar_priors(3, prior_S = matrix(1, 3, 3)), "positive definite")
  expect_error(gvar_priors(0), "K")
})

test_that("stan_gvar warns about unknown priors", {
  skip_on_cran()
  data(ts_data)
  example_data <- ts_data[1:50, 1:3]

  expect_warning(
    stan_gvar(example_data,
      priors = list(prior_s = diag(3)),
      cov_prior = "IW",
      n_chains = 1,
      iter_sampling = 50,
      iter_warmup = 50,
      refresh = 0
    ),
    "Unknown prior"
  )
})

fake_fit <- function(cnames = c("a", "b", "c")) {
  structure(
    list(stan_fit = NULL, arguments = list(n_t = 40, cnames = cnames)),
    class = c("tsnet_fit", "list")
  )
}

test_that("priors_from_fit returns moments in Stan orientation", {
  set.seed(1)
  n_draws <- 2000
  draws <- array(rnorm(n_draws * 9, mean = 0, sd = 0.01), c(n_draws, 3, 3))
  # effect of variable 2 (lagged) on variable 1
  draws[, 1, 2] <- rnorm(n_draws, mean = 0.5, sd = 0.1)
  local_mocked_bindings(.extract_draws = function(fit, par) draws)

  priors <- priors_from_fit(fake_fit())

  expect_named(priors, c("prior_Beta_loc", "prior_Beta_scale"))
  expect_equal(priors$prior_Beta_loc[1, 2], 0.5, tolerance = 0.02)
  expect_equal(priors$prior_Beta_loc[2, 1], 0, tolerance = 0.01)
  expect_equal(priors$prior_Beta_scale[1, 2], 0.1, tolerance = 0.1)
  # sds of 0.01 are floored
  expect_equal(priors$prior_Beta_scale[2, 1], 0.05)
  expect_equal(rownames(priors$prior_Beta_loc), c("a", "b", "c"))

  inflated <- priors_from_fit(fake_fit(), inflate = 3, min_scale = 0.001)
  expect_equal(
    inflated$prior_Beta_scale[1, 2],
    3 * priors$prior_Beta_scale[1, 2]
  )
  expect_lt(inflated$prior_Beta_scale[2, 1], 0.05)

  # usable as input for gvar_priors
  expect_no_error(do.call(gvar_priors, c(list(K = 3), priors)))
})

test_that("priors_from_fit validates its input", {
  expect_error(priors_from_fit(list()), "tsnet_fit")
  expect_error(priors_from_fit(fake_fit(), inflate = 0), "inflate")
  expect_error(priors_from_fit(fake_fit(), min_scale = -1), "min_scale")
})
