#' Compute Effective Sample Sizes for MCMC Samples of BGGM GVAR models
#'
#' This function computes the effective sample sizes (ESS) of MCMC samples of VAR and partial correlation coefficients (pcor) based on the provided MCMC fit object.
#'
#' @param fitobj A list containing a BGGM fit object.
#' @param burnin An integer indicating the number of burn-in iterations to discard. Default is 50.
#'
#' @return A list with two elements: ess_beta and ess_pcor. ess_beta contains the ESS of MCMC samples of VAR, and ess_pcor contains the ESS of MCMC samples of partial correlation coefficients.
#'

#'
#' @import coda
#' @export

ess_gvar <- function(fitobj,
                    burnin = 50){
  # Input Information
  it <- fitobj$iter
  p <- fitobj$p

  ## Get samples
  beta <- fitobj$fit$beta[,,(burnin+1):(iterations+burnin)]
  pcor <- fitobj$fit$pcors[,,(burnin+1):(iterations+burnin)]

  # Transform to mcmc objects
  mcmc_beta <- coda::as.mcmc(t(matrix(beta, p*p, iterations)))
  mcmc_pcor <- coda::as.mcmc(t(matrix(pcor, p*p, iterations)))

  # correct variable names
  # column after column
  cnames <- colnames(fitobj$Y)
  cnames_lag <- paste0(colnames(fitobj$Y), ".l1")

  beta_names <- c(sapply(cnames, paste, cnames_lag, sep = "--"))
  pcor_names <- c(sapply(cnames, paste, cnames, sep = "--"))

  ## Calculate ESS
  ess_beta <- coda::effectiveSize(mcmc_beta)
  ess_pcor <- coda::effectiveSize(mcmc_pcor)

  names(ess_beta) <- beta_names
  names(ess_pcor) <- pcor_names

  ## Return
  l_out <- list(
    ess_beta = ess_beta,
    ess_pcor = ess_pcor
  )
  return(l_out)

}
