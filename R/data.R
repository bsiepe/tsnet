#' Simulated Time Series Dataset
#'
#' This dataset contains a simulated time series dataset for two individuals
#' generated using the graphicalVAR package.
#'
#' @format ## `ts_data`
#' A data frame with 500 rows and 7 columns.
#' @source Simulated using graphicalVAR::graphicalVARsim function.
#'
#' @usage data(ts_data)
#'
#' @details
#' The dataset consists of 250 observations each of 6 variables for two individuals.
#' The variables V1-V6 represent simulated time series data generated using the graphicalVARsim function from the graphicalVAR package.
#' The 'id' column contains a character string as identifier of the two individuals.
#' The data have been standardized using the `scale` function to have zero mean and unit variance.
#'
#' @examples
#' data(ts_data)
#' head(ts_data)
#'
#' @keywords dataset
"ts_data"
