% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_samples.R
\name{stan_fit_convert}
\alias{stan_fit_convert}
\title{Convert Stan Fit to Array of Samples}
\usage{
stan_fit_convert(stan_fit, return_params = c("beta", "sigma", "pcor"))
}
\arguments{
\item{stan_fit}{A Stan fit object obtained from rstan or a tsnet_fit object
from \code{\link{stan_gvar}}.}

\item{return_params}{A character vector specifying which parameters to
return. Options are "beta" (temporal network), "sigma" (innovation
covariance), and "pcor" (partial correlations). Default is
\code{c("beta","sigma", "pcor")}.}
}
\value{
A list containing 3D arrays for the selected parameters. Each array
  represents the posterior samples for a parameter, and each slice of the
  array represents a single iteration.
}
\description{
This function converts a Stan fit object into an array of samples for the
temporal coefficients and the innovation covariance or partial correlation
matrices. It supports \code{rstan} as a backend. It can be used to convert models
fit using \code{\link{stan_gvar}} into 3D arrays, which is the standard data structure
used in \code{tsnet}. The function allows to select which parameters should be
returned.
}
\examples{
\donttest{
data(ts_data)
example_data <- ts_data[1:100,1:3]
fit <- stan_gvar(data = example_data,
                 n_chains = 2,
                 n_cores = 1)
samples <- stan_fit_convert(fit, return_params = c("beta", "pcor"))
}

}
