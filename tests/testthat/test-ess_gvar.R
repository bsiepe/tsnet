test_that("ess_gvar calculates effective sample sizes correctly", {
  data("ts_data")
  data1 <- subset(ts_data, id == "ID1")
  data1 <- subset(data1, select = -c(id))
  fitobj <- BGGM::var_estimate(data1)
  burn <- 50
  it <- fitobj$iter
  beta <- fitobj$fit$beta[,,(burn+1):(it+burn)]
  pcor <- fitobj$fit$pcors[,,(burn+1):(it+burn)]

  # Transform to mcmc objects
  expected_ess_beta <- coda::effectiveSize(coda::as.mcmc(t(matrix(beta, fitobj$p * fitobj$p, it))))
  expected_ess_pcor <- coda::effectiveSize(coda::as.mcmc(t(matrix(pcor, fitobj$p * fitobj$p, it))))

  result <- ess_gvar(fitobj, burnin = burn)

  # Check if the calculated ESS values match the expected values
  # ignore names
  expect_equal(result$ess_beta, expected_ess_beta, ignore_attr = TRUE)
  expect_equal(result$ess_pcor, expected_ess_pcor, ignore_attr = TRUE)
})
