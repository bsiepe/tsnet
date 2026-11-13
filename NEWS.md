# tsnet (development version)

* New `simulate_gvar()` simulates data from a GVAR model with known temporal
  coefficients and innovation covariance matrix.
* New `gvar_stationary()` computes the stationary covariance matrix and the
  standardized parameters implied by a GVAR model, which are the appropriate
  reference when checking parameter recovery of a model fitted to
  standardized data.
* `stan_gvar()` gains the argument `prior_only` to sample from the prior
  distribution without conditioning on the data. `print()` indicates
  prior-only fits and `compare_gvar()` does not accept them.
* New `simulate()` method for `tsnet_fit` objects for posterior predictive
  simulation, or prior predictive simulation for prior-only fits.
* New `gvar_priors()` returns the default priors of `stan_gvar()` with
  optional overrides and checks names and dimensions.
* New `priors_from_fit()` turns the posterior of the temporal coefficients of
  a fitted model into priors for a new model.
* `stan_gvar()` now warns if `priors` contains unknown names, which were
  silently ignored before.

# tsnet 0.2.0
* Fixed documentation issues (e.g., wrong display of LaTex equations).
* Fixed test for ggplot object.
* Added experimental forecasting functionality.

# tsnet 0.1.0

* Initial CRAN submission.
