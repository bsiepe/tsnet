# Orientation round trip: simulate_gvar() -> stan_gvar() -> Stan's Beta.
# Run from the package root: Rscript dev/round_trip.R
# Writes dev/output/round_trip.txt.

devtools::load_all(quiet = TRUE)
dir.create("dev/output", showWarnings = FALSE)

Beta <- matrix(c(
  0.3, 0.5, 0.0,
  0.0, 0.2, 0.0,
  0.0, -0.3, 0.4
), nrow = 3, byrow = TRUE)

round_trip <- function(Sigma, cov_prior, seed) {
  set.seed(seed)
  sim_data <- simulate_gvar(Beta, Sigma, n_t = 5000)
  truth <- gvar_stationary(Beta, Sigma)
  fit <- suppressWarnings(stan_gvar(sim_data,
    cov_prior = cov_prior,
    n_chains = 2,
    iter_sampling = 500,
    iter_warmup = 500,
    refresh = 0,
    seed = seed
  ))
  beta_mean <- apply(rstan::extract(fit$stan_fit, "Beta")$Beta, c(2, 3), mean)
  rho_mean <- apply(rstan::extract(fit$stan_fit, "Rho")$Rho, c(2, 3), mean)
  list(
    truth_std = truth$Beta_std,
    beta_mean = beta_mean,
    max_abs_error = max(abs(beta_mean - truth$Beta_std)),
    max_abs_error_transposed = max(abs(beta_mean - t(truth$Beta_std))),
    max_abs_error_unstandardized = max(abs(beta_mean - Beta)),
    max_abs_error_pcor = max(abs(rho_mean - truth$pcor)),
    # stan_fit_convert() is documented to return the transpose
    convert_is_transpose = isTRUE(all.equal(
      unname(stan_fit_convert(fit, return_params = "beta")$beta_mu),
      unname(t(beta_mean))
    ))
  )
}

Sigma_unequal <- diag(c(1, 4, 0.25))
Sigma_unequal[1, 2] <- Sigma_unequal[2, 1] <- 0.8

res <- list(
  "equal innovation variances, LKJ" = round_trip(diag(3), "LKJ", 1),
  "equal innovation variances, IW" = round_trip(diag(3), "IW", 2),
  "unequal innovation variances, IW" = round_trip(Sigma_unequal, "IW", 3)
)

options(digits = 3)
sink("dev/output/round_trip.txt")
cat("tsnet orientation round trip,", format(Sys.time()), "\n")
cat("n_t = 5000, K = 3, 2 chains x 500 draws\n\ntrue Beta (row = outcome):\n")
print(Beta)
for (nm in names(res)) {
  cat("\n==", nm, "==\n")
  cat("standardized truth:\n")
  print(res[[nm]]$truth_std)
  cat("posterior mean of Stan's Beta:\n")
  print(res[[nm]]$beta_mean)
  print(unlist(res[[nm]][-(1:2)]))
}
sink()
