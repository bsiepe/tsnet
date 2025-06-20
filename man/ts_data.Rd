% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ts_data}
\alias{ts_data}
\title{Simulated Time Series Dataset}
\format{
## `ts_data` A data frame with 500 rows and 7 columns.
\describe{
  \item{id}{A character string identifier for the individual. There are two unique ids, representing two individuals.}
  \item{V1-V6}{These columns represent six different variables in the time series data.}
  }
}
\source{
Simulated using the \code{\link[graphicalVAR]{graphicalVARsim}} function in the \code{graphicalVAR} package.
}
\usage{
data(ts_data)
}
\description{
This dataset contains a simulated time series dataset for two individuals
generated using the \code{graphicalVAR} package. The dataset is useful for testing
and demonstrating the functionality of the package.
}
\details{
The dataset consists of 250 observations each of 6 variables for two individuals.
The variables V1-V6 represent simulated time series data generated using the graphicalVARsim function from the graphicalVAR package.
The 'id' column contains a character string as identifier of the two individuals.
The data have been standardized  to have zero mean and unit variance.
}
\keyword{dataset}
