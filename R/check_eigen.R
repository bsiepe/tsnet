#' Check Eigenvalues of Bayesian GVAR object
#'
#' This function checks the eigenvalues of the Beta matrix (containing the
#' temporal coefficients) to assure that the model is stationary. It uses the
#' same check as the `graphicalVAR` package. The function calculates the
#' eigenvalues of the Beta matrix and checks if the sum of the squares of the
#' real and imaginary parts of the eigenvalues is less than 1. If it is, the VAR
#' model is considered stable.
#'
#' @param fitobj A fitted Bayesian GVAR object. This can be a tsnet_fit object
#'   (obtained from [stan_gvar()]), a BGGM object (obtained from
#'   [BGGM::var_estimate()]), or extracted posterior samples (obtained from
#'   [stan_fit_convert()).
#' @param verbose Logical. If TRUE, a verbal summary of the results is printed.
#'   Default is TRUE.
#' @examples
#'  data(fit_data)
#'  fitobj <- fit_data[[1]]
#'  result <- check_eigen(fitobj)
#'
#' @return A list containing the eigenvalues and a verbal summary of the
#'   results.
#' @export
check_eigen <- function(fitobj,
                        verbose = TRUE){

  # Input check
  if(!(inherits(fitobj, "var_estimate") ||
       inherits(fitobj, "tsnet_fit") ||
       inherits(fitobj, "tsnet_samples"))) {
    stop("Error: 'fitobj' must be either a 'var_estimate', 'tsnet_fit', or 'tsnet_samples' object.")
  }

  if(inherits(fitobj, "tsnet_fit")) {
    fitobj <- stan_fit_convert(fitobj,
                                  return_params = c("beta", "pcor"))
  }

  # Extract elements
  beta_mu <- fitobj$beta_mu
  # Calculate eigenvalues
  eigen_beta_mu <- eigen(beta_mu)$values

  # check if the sum of the squares of the real and imaginary parts of the eigenvalues of the coefficient matrix is less than 1
  # If it is, it implies that the VAR model is stable
  check_beta_mu <- all(Re(eigen_beta_mu)^2 + Im(eigen_beta_mu)^2 <1)
  if(isTRUE(check_beta_mu)){
    response <- "stable."
  }
  else {
     response <- "not stable. It is recommended to refit the model and to inspect your data."
  }

  # Return verbal summary
  if(verbose) {
    message(paste0("The VAR coefficient matrix of the input model is ", response))
  }

  # Return results
  list(eigenvalues = list(eigen_beta_mu = eigen_beta_mu))
}

