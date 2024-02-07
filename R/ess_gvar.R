#' Compute Effective Sample Sizes for MCMC Samples of BGGM GVAR models
#'
#' This function computes the effective sample size (ESS) of MCMC samples of
#' temporal and contemporaneous network parameters based on the provided
#' [BGGM::var_estimate()] fit object. It uses the default functionality of the
#' `coda` package to compute ESS. Effective sample size estimates for
#' [stan_gvar()] models can be found in their model output itself. Currently
#' used for internal development only.
#'
#' @param fitobj A [BGGM::var_estimate()] fit object.
#' @param burnin An integer indicating the number of burn-in iterations to
#'   discard. Default is 0.
#'
#' @return A list with two elements: ess_beta and ess_pcor. ess_beta contains
#'   the ESS of MCMC samples of the temporal coefficients, and ess_pcor contains
#'   the ESS of MCMC samples of partial correlation coefficients.
#'
#' @examples
#' \dontrun{
#' # Load data
#' data(ts_data)
#' example_data <- ts_data[1:100,]
#'
#' # Estimate a GVAR model
#' fit <- BGGM::var_estimate(example_data)
#' ess <- ess_gvar(fit)
#' }
#'
#' @importFrom coda as.mcmc effectiveSize
#' @noRd

ess_gvar <- function(fitobj,
                     burnin = 0) {


  # Input Information
  it <- fitobj$iter
  p <- fitobj$p

  ## Get samples
  beta <- fitobj$fit$beta[, , (burnin + 1):(it + burnin)]
  pcor <- fitobj$fit$pcors[, , (burnin + 1):(it + burnin)]

  # Transform to mcmc objects
  mcmc_beta <- coda::as.mcmc(t(matrix(beta, p * p, it)))
  mcmc_pcor <- coda::as.mcmc(t(matrix(pcor, p * p, it)))

  # correct variable names
  # column after column
  cnames <- colnames(fitobj$Y)
  cnames_lag <- paste0(colnames(fitobj$Y), ".l1")

  beta_names <- as.vector(vapply(cnames,
                                  function(x) paste(x, cnames_lag, sep = "--"),
                                  character(length(cnames))))
  pcor_names <- as.vector(vapply(cnames,
                                 function(x) paste(x, cnames, sep = "--"),
                                 character(length(cnames))))

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
