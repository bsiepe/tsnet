# test_that("ess_gvar calculates effective sample sizes correctly", {
#   data("ts_data")
#   example_data <- ts_data[1:100,1:3]
#   fit <- stan_gvar(example_data, n_chains = 1, iter_sampling = 500)
#   it <- 500
#   p <- 3
#   samples <- stan_fit_convert(fit)
#
#
#   # Transform to mcmc objects
#   expected_ess_beta <- coda::effectiveSize(coda::as.mcmc(t(matrix(samples$fit$beta, p*p, it))))
#   expected_ess_pcor <- coda::effectiveSize(coda::as.mcmc(t(matrix(samples$fit$pcor, p*p, it))))
#
#   result <- ess_gvar(fit, burnin = 0)
#
#   # Check if the calculated ESS values match the expected values
#   # ignore names
#   expect_equal(result$ess_beta, expected_ess_beta, ignore_attr = TRUE)
#   expect_equal(result$ess_pcor, expected_ess_pcor, ignore_attr = TRUE)
# })
