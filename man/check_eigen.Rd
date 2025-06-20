% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_eigen.R
\name{check_eigen}
\alias{check_eigen}
\title{Check Eigenvalues of Bayesian GVAR object}
\usage{
check_eigen(fitobj, verbose = TRUE)
}
\arguments{
\item{fitobj}{A fitted Bayesian GVAR object. This can be a tsnet_fit object
(obtained from \code{\link{stan_gvar}}, a \code{BGGM} object (obtained from
\code{\link[BGGM]{var_estimate}}), or extracted posterior samples
(obtained from \code{\link{stan_fit_convert}}.}

\item{verbose}{Logical. If TRUE, a verbal summary of the results is printed.
Default is \code{TRUE}.}
}
\value{
A list containing the eigenvalues and a verbal summary of the
  results.
}
\description{
This function checks the eigenvalues of the Beta matrix (containing the
temporal coefficients) to assure that the model is stationary. It uses the
same check as the \code{graphicalVAR} package. The function calculates the
eigenvalues of the Beta matrix and checks if the sum of the squares of the
real and imaginary parts of the eigenvalues is less than 1. If it is, the VAR
model is considered stable.
}
\examples{
 data(fit_data)
 fitobj <- fit_data[[1]]
 result <- check_eigen(fitobj)

}
