# tsnet (development version)

* The default of `prior_Rho_scale` for `cov_prior = "LKJ"` changed from
  `sqrt(0.5)` to `2`. `prior_Rho_scale` is the precision of a beta prior on
  the rescaled partial correlations. The new default is the uniform prior
  that the documentation described, whereas the previous default was
  U-shaped and favored partial correlations close to -1 and 1. Results of
  models with `cov_prior = "LKJ"` and default priors therefore change; pass
  `priors = list(prior_Rho_scale = matrix(sqrt(0.5), K, K))` to reproduce
  results of version 0.2.0. Models with `cov_prior = "IW"` are not affected.
* The documentation of the priors in `stan_gvar()` was corrected: the
  default of `cov_prior` is `"IW"`, the inverse Wishart prior is placed on
  the precision matrix, and `prior_Rho_marginal` is the prior variance (not
  the standard deviation) of the partial correlations.
* New `simulate_gvar()` simulates data from a GVAR model with known temporal
  coefficients and innovation covariance matrix.
* New `gvar_stationary()` computes the stationary covariance matrix and the
  standardized parameters implied by a GVAR model, which are the appropriate
  reference when checking parameter recovery of a model fitted to
  standardized data.
* `stan_gvar()` gains the argument `prior_only` to sample from the prior
  distribution without conditioning on the data. `print()` indicates
  prior-only fits and `compare_gvar()` does not accept them. The values of
  the data are not used for prior-only fits, so placeholder data can be
  passed.
  Prior-only fits with `cov_prior = "LKJ"` give a warning if the beta prior
  on a partial correlation is unbounded.
* New `simulate()` method for `tsnet_fit` objects for posterior predictive
  simulation, or prior predictive simulation for prior-only fits. By
  default, non-stationary draws are not used, and a message reports how
  many were dropped.
* New `gvar_priors()` returns the default priors of `stan_gvar()` with
  optional overrides and checks names and dimensions.
* New `priors_from_fit()` turns the posterior of the temporal coefficients of
  a fitted model into priors for a new model.
* New `prior_posterior_summary()` reports prior and posterior means and
  standard deviations as well as the posterior contraction for each edge.
* `stan_gvar()` now warns if `priors` contains unknown names, which were
  silently ignored before.
* `stan_gvar()` now stores the sample means and standard deviations used to
  standardize the data in `arguments$center` and `arguments$scale`.

# tsnet 0.2.0
* Fixed documentation issues (e.g., wrong display of LaTex equations).
* Fixed test for ggplot object.
* Added experimental forecasting functionality.

# tsnet 0.1.0

* Initial CRAN submission.
