#' Check Eigenvalues of BGGM Object
#'
#' This function checks the eigenvalues of Beta to assure that the model is stationary. It uses the same check as the `graphicalVAR` package. 
#'
#' @param fitobj A BGGM object.
#' @return A list containing the eigenvalues and a verbal summary of the results.
check_eigen <- function(fitobj){
  # Extract elements
  beta_mu <- fitobj$beta_mu
  # Calculate eigenvalues
  eigen_beta_mu <- eigen(beta_mu)$values
  
  # check if the sum of the squares of the real and imaginary parts of the eigenvalues of the coefficient matrix is less than 1 
  # If it is, it implies that the VAR model is stable
  check_beta_mu <- all(Re(eigen_beta)^2 + Im(eigen_beta)^2 <1)
  if(isTRUE(check_beta_mu)){
    response <- "stable."
  }
  else {
     response <- "not stable. It is recommended to refit the model and to inspect your data."
  }
  
  # Return verbal summary
  print(paste0("The VAR coefficient matrix of the input model is ", response))

  # Return results
  list(eigenvalues = list(eigen_beta_mu = eigen_beta_mu))
}
