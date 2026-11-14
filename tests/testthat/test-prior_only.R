test_that("prior_only is validated before sampling", {
  data(ts_data)
  example_data <- ts_data[1:50, 1:3]

  expect_error(
    stan_gvar(example_data, prior_only = TRUE, ahead = 1),
    "ahead"
  )
  expect_error(
    stan_gvar(example_data, prior_only = TRUE, method = "variational"),
    "sampling"
  )
})

test_that("unbounded beta priors on partial correlations are detected", {
  expect_false(.unbounded_rho_prior(gvar_priors(3)))
  # default up to version 0.2.0
  expect_true(.unbounded_rho_prior(
    gvar_priors(3, prior_Rho_scale = matrix(sqrt(0.5), 3, 3))
  ))
  # only the upper triangle is used by the model
  scale <- matrix(2, 3, 3)
  scale[2, 1] <- 0.5
  expect_false(.unbounded_rho_prior(gvar_priors(3, prior_Rho_scale = scale)))
  scale[1, 2] <- 0.5
  expect_true(.unbounded_rho_prior(gvar_priors(3, prior_Rho_scale = scale)))
  loc <- matrix(0.5, 3, 3)
  loc[1, 3] <- 0.9
  expect_true(.unbounded_rho_prior(
    gvar_priors(3, prior_Rho_loc = loc, prior_Rho_scale = matrix(5, 3, 3))
  ))
})

test_that("prior-only LKJ fits with unbounded beta priors give a warning", {
  skip_on_cran()
  data(ts_data)

  # sampler warnings are expected here and not of interest
  suppressWarnings(expect_warning(
    stan_gvar(ts_data[1:50, 1:3],
      priors = list(prior_Rho_scale = matrix(sqrt(0.5), 3, 3)),
      cov_prior = "LKJ",
      n_chains = 1,
      iter_sampling = 20,
      iter_warmup = 20,
      refresh = 0,
      prior_only = TRUE
    ),
    "Prior-only sampling"
  ))
})

test_that("prior_only fits sample from the prior", {
  skip_on_cran()
  data(ts_data)
  example_data <- ts_data[1:50, 1:3]
  loc <- matrix(c(
    0.5, -0.4, 0,
    0, 0.2, 0,
    0, 0, 0
  ), nrow = 3, byrow = TRUE)

  for (cov_prior in c("IW", "LKJ")) {
    fit <- suppressWarnings(stan_gvar(example_data,
      priors = list(prior_Beta_loc = loc),
      cov_prior = cov_prior,
      n_chains = 1,
      iter_sampling = 500,
      iter_warmup = 300,
      refresh = 0,
      seed = 1,
      prior_only = TRUE
    ))

    expect_true(fit$arguments$prior_only)
    beta_draws <- rstan::extract(fit$stan_fit, "Beta")$Beta
    # prior sd is 0.5, the data would pull estimates towards much less
    expect_lt(max(abs(apply(beta_draws, c(2, 3), mean) - loc)), 0.12)
    expect_lt(max(abs(apply(beta_draws, c(2, 3), sd) - 0.5)), 0.12)
  }

  expect_output(print(fit), "prior-only fit")
  expect_error(compare_gvar(fit, fit), "prior-only")
  expect_no_error(stan_fit_convert(fit))
  expect_no_error(suppressMessages(check_eigen(fit)))
  expect_no_error(suppressWarnings(posterior_plot(fit)))

  sims <- suppressWarnings(simulate(fit, nsim = 20, seed = 1))
  expect_true(all(attr(sims, "draws")$stationary))
  expect_named(sims[[1]], colnames(example_data))
})

test_that("regular fits are not flagged as prior-only", {
  skip_on_cran()
  data(ts_data)
  example_data <- ts_data[1:50, 1:3]
  fit <- suppressWarnings(stan_gvar(example_data,
    cov_prior = "IW",
    n_chains = 1,
    iter_sampling = 100,
    iter_warmup = 100,
    refresh = 0
  ))

  expect_false(fit$arguments$prior_only)
  expect_false(any(grepl("prior-only", capture.output(print(fit)))))
})

test_that("simulated data are recovered in the correct orientation", {
  skip_on_cran()
  Beta <- matrix(c(
    0.3, 0.5, 0.0,
    0.0, 0.2, 0.0,
    0.0, -0.3, 0.4
  ), nrow = 3, byrow = TRUE)
  Sigma <- diag(3)
  set.seed(1)
  sim_data <- simulate_gvar(Beta, Sigma, n_t = 1000)
  truth <- gvar_stationary(Beta, Sigma)$Beta_std

  fit <- suppressWarnings(stan_gvar(sim_data,
    cov_prior = "LKJ",
    n_chains = 1,
    iter_sampling = 300,
    iter_warmup = 300,
    refresh = 0,
    seed = 1
  ))
  beta_mean <- apply(rstan::extract(fit$stan_fit, "Beta")$Beta, c(2, 3), mean)

  expect_lt(max(abs(beta_mean - truth)), 0.12)
  expect_gt(max(abs(beta_mean - t(truth))), 0.3)
  # stan_fit_convert() returns the transposed matrix
  expect_equal(
    unname(stan_fit_convert(fit, return_params = "beta")$beta_mu),
    unname(t(beta_mean))
  )
})
