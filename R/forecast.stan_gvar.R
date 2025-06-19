#' Forecast Using Fitted GVAR Models
#'
#' @description
#' This function generates forecasts using one or more fitted GVAR models.
#' It supports forecasting with model combinations, optionally weighted by model performance,
#' and handles scenarios with or without new data for prediction.
#'
#' @param fitobj A \code{tsnet_fit} object or a \code{tsnet_list} object containing one or more fitted GVAR models.
#' @param ahead An integer specifying the forecast horizon. Default is \code{1}.
#' @param new_data Optional. A matrix or array of new data for forecasting. If not provided, forecasts are generated for an empty array.
#' @param interval A numeric vector specifying the prediction interval (e.g., \code{c(0.025, 0.975)}). Values must lie between 0 and 1.
#' @param weights An object of class \code{model_weights}, typically obtained from \code{\link{model_weights}}.
#' If not provided, the function computes model weights internally.
#' @param start Integer indicating the starting point for cross-validation. Inherited from \code{\link{lfo.stan_gvar}}. Default is \code{1}.
#' @param method A character string specifying the method for model combination. Default is \code{"stacking"}.
#' @param inc_samples Logical. If \code{TRUE}, includes all posterior samples in the output. Default is \code{FALSE}.
#' @param ... Additional arguments (currently unused).
#'
#' @details
#' The function processes the input models (either a single \code{tsnet_fit} or a list of them),
#' validates that all models are based on the same time series, and prepares data for forecasting.
#' If no model weights are provided, the function computes them using \code{\link{model_weights}}.
#' Forecasts are computed for the specified horizon (\code{ahead}) and optionally weighted by model performance.
#'
#' The forecasts include posterior predictive intervals based on the specified \code{interval} argument.
#' If \code{new_data} is not provided, forecasts assume no additional data and generate predictions accordingly.
#'
#' @return A \code{tsnet_forecast} object with the following components:
#' \describe{
#'   \item{forecast}{A list containing the forecasted means and optionally the full posterior samples.}
#'   \item{mean}{The mean forecast values.}
#'   \item{samples}{(Optional) A list of forecasted posterior samples if \code{inc_samples} is \code{TRUE}.}
#' }
#'
#'
#' @examples
#' \dontrun{
#' # Assuming `fit1` and `fit2` are tsnet_fit objects
#' fits <- tsnet_list(fit1, fit2)
#' forecasts <- forecast.stan_gvar(fits, ahead = 5, method = "stacking")
#' print(forecasts)
#' }
#' @keywords internal
#' @noRd


forecast.stan_gvar <- function(fitobj,
                               ahead = 1,
                               new_data = NULL,
                               interval = c(0.025, 0.975),
                               weights = NULL,
                               start = 1,
                               method = "stacking",
                               inc_samples = FALSE,
                               ...) {

  ## Input checks
  # Input can either be a single model or a list of models
  # check if the input is a tsnet_fit or a tsnet_list
  if (!inherits(fitobj, "tsnet_fit") &
      !inherits(fitobj, "tsnet_list")) {
    stop("Input must be a tsnet_fit or a tsnet_list object")
  }

  # Check if the interval is between 0 and 1
  if (any(interval < 0) | any(interval > 1)) {
    stop("Interval must be between 0 and 1")
  }
  # Check if weights are specified correctly
  if(!is.null(weights) && !inherits(weights, "model_weights")){
    stop("Model weights must be obtained from the model_weights function.")
  }


  ## Prepare data
  # number of models
  n_mods <- 1
  if (inherits(fitobj, "tsnet_list")) {
    n_mods <- length(fitobj)
  }
  else {
    fitobj <- tsnet_list(fitobj)
  }

  # check if they all used the same time series
  col_names <- fitobj[[1]]$arguments$cnames
  colnames_ident <- all(sapply(fitobj, function(x){
    all(x$names == col_names)
  }))
  if(!colnames_ident){
    stop("All models must use the same time series.")
  }

  # newdata available?
  if (is.null(new_data)) {
    newdata <- array(0, dim = c(ahead, fitobj[[1]]$arguments$n_t))
    compute_log_lik <- 0
  } else {
    compute_log_lik <- 1
  }

  # model weights available?
  # Case 1: No model weights provided
  # Case 2: model weights from a model_weight object
  # Case 3: No model weights requested
  if(n_mods > 1 & is.null( weights ) ) {
    mw <- model_weights(fitobjs = fitobj, start = start)
    weights <- mw$weights []
  } else if( n_mods > 1 & !is.null( weights ) ) {
    weights <- weights$weights[]
  } else if( n_mods == 1 ) {
    weights <- 1
  }


  ## Compute Forecasts
  # Looping over all specified models
  forecast_obj <- lapply(fitobj, function(m){

    args <- m$arguments$fn_args
    prior_delta <- (1 / m$arguments$priors$prior_Rho_marginal) - 1

    standat <- list(Y = as.matrix(m$arguments$data),
                    K = m$arguments$p,
                    "T" = m$arguments$n_t,
                    beep = eval(args$beep),
                    prior_Beta_loc = m$arguments$priors$prior_Beta_loc,
                    prior_Beta_scale = m$arguments$priors$prior_Beta_scale,
                    prior_Rho_loc = m$arguments$priors$prior_Rho_loc,
                    prior_Rho_scale = m$arguments$priors$prior_Rho_scale,
                    prior_Eta = m$arguments$priors$prior_Eta,
                    prior_S = m$arguments$priors$prior_S,
                    prior_delta = prior_delta,
                    ahead = ahead,
                    Y_future = new_data,
                    compute_log_lik = compute_log_lik)

  # Obtain model name for sampling
  gqs_model <- .get_gvar_name(m)

  # backcast <- max(m$mgarchP, m$mgarchQ)
  # this is just lag-1 here
  backcast <- 1

  # TODO think really good about the indexing here
  nt <- m$arguments$n_t
  cast_start <- (nt - backcast + 1)
  forecast_start <- (nt + 1)
  forecast_end <- (nt + ahead)

  forecasted <- rstan::gqs(stanmodels[[gqs_model]],
                           draws = as.matrix(m$stan_fit),
                           data = standat)
  return(forecasted)
  })

  ## Extract forecasts
  # extract forecasted values with weighting
  # f_mean <-
  # restructure to correct format




  ## Output
  # Return a tsnet_forecast object
  out <- list()
  out$forecast <- list()
  out$forecast$forecast_obj <- forecast_obj
  # out$forecast$weighted_forecasts <- .weighted_samples(forecast_obj,
  #                                                      params = "Y_forecast",
  #                                                      weights = weights)


  # should all samples be included?
  # if(inc_samples){
    # f_samples <- .weighted_samples(forecast_obj, weights = weights, method = method)
    # out$forecast$samples$mean <- f_samples
  # }
  if(compute_log_lik == 1 ) {
    log_lik <- lapply(forecast_obj, function(x) {
      rstan::extract(x, pars = "log_lik")$log_lik
    })
    out$forecast$log_lik <- log_lik
  }


  # Store arguments
  mc <- match.call()
  fl <- formals()
  missing_args <- setdiff(names(fl), names(mc))
  missing_args <- Map(function(arg, default)
    if (!is.null(default)) mc[[arg]] <- default, missing_args, fl[missing_args])
  all_args <- c(as.list(mc), missing_args)

  args <- list(
    ahead = ahead,
    new_data = new_data,
    interval = interval,
    weights = weights,
    start = start,
    method = method,
    inc_samples = inc_samples,
    fn_args = all_args
  )

  out$arguments <- args

  class(out) <- c("tsnet_forecast", class(out))
  return(out)

}




# Helpers

#' Create a List of \code{tsnet_fit} Objects
#'
#' @description
#' This internal function creates a list of \code{tsnet_fit} objects, which are the
#' results of fitting a time series network model using the\code{\link{stan_gvar}} function.
#' The list is assigned the class \code{tsnet_list}, making it compatible with other
#' functions that expect such a list.
#'
#' @param ... Objects of class \code{tsnet_fit}, which are the results of fitting
#' time series network models. These objects are passed as arguments to the function.
#'
#' @details
#' The function takes any number of \code{tsnet_fit} objects (from the \code{\link{stan_gvar}}
#' model fitting process), combines them into a list, and assigns the class \code{tsnet_list}.
#' This allows the list to be used as input in other functions that expect a \code{tsnet_list}
#' object, providing compatibility with workflows that involve multiple time series models.
#'
#' @return A list of \code{tsnet_fit} objects, with the class \code{tsnet_list}.
#'
#' @seealso \code{\link{stan_gvar}}, \code{\link{forecast.stan_gvar}}, \code{\link{model_weights}}
#'
#' @keywords internal
#' @noRd

tsnet_list <- function(...) {
  out <- list(...)
  class(out) <- c("tsnet_list", class(out))
  return(out)
}

# Obtain model name from stan_gvar output

#' Get gvar name
#'
#' Obtain the name of the \code{\link{stan_gvar}} model from a fitte model output.
#'
#' @param fitobj A \code{\link{stan_gvar}} object
#'
#' @return character name of the \code{\link{stan_gvar}} model
#'
#' @keywords internal
#' @noRd
.get_gvar_name <- function(fitobj){
  if(fitobj$arguments$fn_args$cov_prior == "IW"){
    if(isTRUE(fitobj$arguments$fn_args$rmv_overnight)){
      gqs_model <- "VAR_wishart_beep"
    }
    else{
      gqs_model <- "VAR_wishart"
    }
  }
  if(fitobj$arguments$fn_args$cov_prior == "LKJ"){
    if(isTRUE(fitobj$arguments$fn_args$rmv_overnight)){
      gqs_model <- "VAR_LKJ_beep"
    }
    else{
      gqs_model <- "VAR_LKJ"
    }
  }
  return(gqs_model)
}




#'
#' @description
#' This internal function applies weights to the posterior samples of different models
#' and sums them together. The function is used to combine the posterior samples
#' from multiple models based on their weights.
#'
#' @param fitobj A list of fitted model objects from which to extract posterior samples.
#' @param params A character vector specifying the parameters to extract from the models.
#' @param weights A numeric vector of weights to apply to the posterior samples of each model.
#'
#' @details
#' The function first extracts the posterior samples for the specified parameters from each model.
#' It then applies the weights to the samples and sums the weighted samples together.
#' The resulting combined samples have the same dimensionality as the original samples.
#'
#' @return A list of combined posterior samples for each specified parameter.
#'
#' @examples
#' \dontrun{
#' # Assuming `fit1` and `fit2` are fitted model objects
#' fitobj <- tsnet_list(fit1, fit2)
#' params <- c("param1", "param2")
#' weights <- c(0.6, 0.4)
#' combined_samples <- .weighted_samples(fitobj, params, weights)
#' }
#' @keywords internal
#' @noRd
.weighted_samples <- function(fitobj, params, weights) {
  samps <- lapply(fitobj, rstan::extract, pars = params)
  for(i in seq_len(length(samps))) { # each model
    samps[[i]] <- lapply(samps[[i]], function(p) { # each parameter
      p * weights[i]
    })
  }
  samps_comb <- lapply(params, function(p) { # each parameter
    Reduce("+", lapply(samps, function(m) {m[[p]]}))
  })
  names(samps_comb) <- params
  return(samps_comb)
}
