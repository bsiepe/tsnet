#' Fit Bayesian Generalized Vector Autoregressive (gVAR) Models with Stan
#'
#' This function fits a Bayesian gVAR model to the provided data using Stan.
#'
#' @param data A data frame or matrix containing the time series data.
#' @param beep A vector of beeps with length of `nrow(data)`. The beep indicator can be used to remove overnight effects from the last beep of a day to the first beep of the next day.
#' @param priors A list of prior distributions for the model parameters.
#' @param method A string indicating the method to use for fitting the model. Options are "sampling" or "variational".
#' @param cov_prior A string indicating the prior distribution to use for the covariance matrix. Options are "LKJ" or "IW" (Inverse-Wishart).
#' @param rmv_overnight A logical indicating whether to remove overnight effects. Default is `FALSE`.
#' @param iter_sampling An integer specifying the number of iterations for the sampling method. Default is 500.
#' @param iter_warmup An integer specifying the number of warmup iterations for the sampling method. Default is 500.
#' @param n_chains An integer specifying the number of chains for the sampling method. Default is 4.
#' @param n_cores An integer specifying the number of cores to use for parallel computation. Default is 4.
#' @param center_only A logical indicating whether to only center (and not scale) the data. Default is `FALSE`.
#' @param ... Additional arguments passed to the `rstan::sampling` or `rstan::vb` function.
#'
#' @return A `stanfit` object representing the fitted model.
#'
#' @examples
#' \dontrun{
#' data <- matrix(rnorm(1000), ncol = 10)
#' fit <- stan_gvar(data, method = "sampling", cov_prior = "LKJ")
#' }
#'
#' @export

stan_gvar <-
  function(data,
           beep = NULL,
           priors = NULL,
           method = "sampling",
           cov_prior = "LKJ",
           rmv_overnight = FALSE,
           iter_sampling = 500,
           iter_warmup = 500,
           n_chains = 4,
           n_cores = 4,
           center_only = FALSE,
           ...) {
    if(isTRUE(center_only)){
      Y <- apply(data, MARGIN = 2, scale, center = TRUE, scale = FALSE)
    } else{
      Y <- apply(data, MARGIN = 2, scale, center = TRUE, scale = TRUE)
    }

    K <- ncol(data)
    n_t <- nrow(data)

    if(is.null(beep)){
      beep <- seq(1, n_t)
    }


    # Specify Priors
    if(is.null(priors[["prior_Beta_loc"]])){
        prior_Beta_loc <- matrix(rep(0, K*K), nrow = K, ncol = K)
    } else {
        prior_Beta_loc <- priors[["prior_Beta_loc"]]
    }

    if(is.null(priors[["prior_Beta_scale"]])){
        prior_Beta_scale <- matrix(rep(.5, K*K), nrow = K, ncol = K)
    } else {
        prior_Beta_scale <- priors[["prior_Beta_scale"]]
    }

    if(is.null(priors[["prior_Rho_loc"]])){
        prior_Rho_loc <- matrix(rep(.5, K*K), nrow = K, ncol = K)
    } else {
        prior_Rho_loc <- priors[["prior_Rho_loc"]]
    }

    if(is.null(priors[["prior_Rho_scale"]])){
        prior_Rho_scale <- matrix(rep(sqrt(.5), K*K), nrow = K, ncol = K)
    } else {
        prior_Rho_scale <- priors[["prior_Rho_scale"]]
    }

    if(is.null(priors[["prior_Rho_marginal"]])){
        prior_Rho_marginal <- 0.25
    } else {
        prior_Rho_marginal <- priors[["prior_Rho_marginal"]]
    }

    if(is.null(priors[["prior_Eta"]])){
      prior_Eta <- 1
    } else {
      prior_Eta <- priors[["prior_Eta"]]
    }

    # Convert SD to delta: SD = 1/(delta+1), delta = (1 / SD) - 1
    prior_delta <- (1 / prior_Rho_marginal) - 1



    # Choose model to fit
    if (cov_prior == "LKJ") {
      # Stan Data
      stan_data <- list(
        K = K,
        "T" = n_t,
        Y = as.matrix(Y),
        beep = beep,
        prior_Rho_loc = prior_Rho_loc,
        prior_Rho_scale = prior_Rho_scale,
        prior_Beta_loc = prior_Beta_loc,
        prior_Beta_scale = prior_Beta_scale,
        prior_Eta = prior_Eta
      )

      if (isTRUE(rmv_overnight)) {
        # remove overnight effects
        model_name <- "VAR_LKJ_beep"
      } else{
        # standard model
        model_name <- "VAR_LKJ"
      }
    }
    if (cov_prior == "IW") {
      # Stan Data
      stan_data <- list(
        K = K,
        "T" = n_t,
        Y = as.matrix(Y),
        beep = beep,
        prior_Rho_loc = prior_Rho_loc,
        prior_Rho_scale = prior_Rho_scale,
        prior_Beta_loc = prior_Beta_loc,
        prior_Beta_scale = prior_Beta_scale,
        prior_delta = prior_delta
      )

      if (isTRUE(rmv_overnight)) {
        # remove overnight effects
        model_name <- "VAR_wishart_beep"
      } else{
        # standard model
        model_name <- "VAR_wishart"
      }
    }

      if (method == "sampling") {
        # Run sampler
        stan_fit <- rstan::sampling(
          object = stanmodels[[model_name]],
          data = stan_data,
          chains = n_chains,
          cores = n_cores,
          iter = iter_sampling + iter_warmup,
          warmup = iter_warmup,
          refresh = 500,
          thin = 1,
          init = .1,
          control = list(adapt_delta = .8),
          ...
        )
      }
      if (method == "variational") {
        stan_fit <- rstan::vb(
          object = stanmodels[[model_name]],
          data = stan_data,
          init = .1,
          tol_rel_obj = .001,
          output_samples = iter_sampling * n_chains,
          ...
        )
      }

    return(stan_fit)
}
