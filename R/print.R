#' Print method for compare_gvar objects
#'
#' This function prints a summary of the Norm-Based Comparison Test for a [compare_gvar()] object.
#'
#' @param x A test object obtained from [compare_gvar()]
#' @param ... Additional arguments to be passed to the print method. (currently not used)
#'
#' @return Prints a summary of the Norm-Based Comparison Test to the console
#'
#' @details This function prints a summary of the Norm-Based Comparison Test for a [compare_gvar()] object.
# It displays the general summary and model-specific results, including the number of significant comparisons
#' in the temporal and contemporaneous networks, as well as the number of reference distances that were larger
#' than the empirical distance for each network.
#'
#' @examples
#' # Load example fits
#' data(fit_data)
#'
#' # Perform test
#' test_res <- compare_gvar(fit_data[[1]], fit_data[[2]], n_draws = 100)
#'
#' # Print results
#' print(test_res)
#'
#' @export
print.compare_gvar <- function(x,
                                 ...){

  if(!inherits(x, "compare_gvar")){
    stop("This function only works with a result of the compare_gvar function.")
  }
  cat("### Summary of the Norm-Based Comparison Test ###")
  cat("\n")
  cat("\n#--- General Summary ---#")
  cat(
    "\nIn the temporal network", x$sig_beta, "of the 2 comparisons were significant."
  )
  cat(
    "\nIn the contemporaneous network", x$sig_pcor, "of the 2 comparisons were significant."
  )
  cat("\n")
  cat(
    "\n#--- Model-specific Results ---#"
  )
  cat(
    "\nFor", x$larger_beta$mod[1], x$larger_beta$sum_larger[1],
    "of the reference distances of the temporal network and",
    x$larger_pcor$sum_larger[1],
    "of the reference distances of the contemporaneous network were larger than the empirical distance."
  )
  cat("\n")
  cat(
    "\nFor", x$larger_beta$mod[2], x$larger_beta$sum_larger[2],
    "of the reference distances of the temporal network and",
    x$larger_pcor$sum_larger[2],
    "of the reference distances of the contemporaneous network were larger than the empirical distance."
  )

}


#' Print method for tsnet_fit objects
#'
#' This method provides a summary of the Bayesian GVAR model fitted with [stan_gvar()].
#' It prints general information about the model, including the estimation method and the number of chains and iterations
#' It also prints the posterior mean of the temporal and contemporaneous coefficients.
#'
#' @param x A tsnet_fit object.
#' @param ... Additional arguments passed to the print method (currently not used).
#'
#' @return Prints a summary to the console.
#'
#' @examples
#' # Load example data
#' data(ts_data)
#' example_data <- ts_data[1:100,1:3]
#'
#' # Fit the model
#' fit <- stan_gvar(example_data,
#'                  method = "sampling",
#'                  cov_prior = "IW",
#'                  n_chains = 2)
#'
#' print(fit)
#' @export
print.tsnet_fit <- function(x,
                            ...) {
  cat("### Summary of the Bayesian GVAR model ###")
  cat("\n")
  cat("\n#--- General Summary ---#")
  # Summary for MCMC
  if (x$arguments$fn_args$method == "sampling") {
    cat(
      "\nModel was estimated using MCMC with",
      x$arguments$fn_args$n_chains,
      "chains using",
      x$arguments$fn_args$iter_sampling,
      "iterations each.Warmup was set to",
        x$arguments$fn_args$iter_warmup,
        "iterations.")
    cat(
      "\nThe model was estimated using the",
      x$arguments$fn_args$cov_prior,
      "covariance prior."
    )
  }

  # Summary for variational inference
  if (x$arguments$fn_args$method == "variational") {
    cat(
      "\nModel was estimated using variational inference with",
      x$arguments$fn_args$iter_sampling * x$arguments$fn_args$n_chains,
      "iterations."
    )
    cat("\nThe model was estimated using the",
        x$arguments$fn_args$cov_prior,
        "covariance prior.")
  }
  cat("\n")

  post_samps <- stan_fit_convert(x, return_params = c("beta", "pcor"))

  cat("\n#--- Parameter Summary ---#")
  cat("\n")
  cat("The posterior mean of the temporal coefficients is:")
  cat("\n")
  print(post_samps$beta_mu)
  cat("\n")
  cat("Rownames correspond to the independent variable, column names to the dependent variable.")
  cat("\n")
  cat("\n")
  cat("The posterior mean of the contemporaneous coefficients is:")
  cat("\n")
  print(post_samps$pcor_mu)
  cat("\n")

}



