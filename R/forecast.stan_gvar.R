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
  if(!inherits(weights, "model_weights")){
    stop("Model weights must be obtained from the model_weights function.")
  }

  # TODO add more checks


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
  col_names <- object[[1]]$arguments$cnames
  colnames_ident <- all(sapply(fitobj, function(x){
    all(x$names == col_names)
  }))
  if(!colnames_ident){
    stop("All models must use the same time series.")
  }

  # newdata available?
  # TODO not implemented yet in Stan
  if (is.null(newdata)) {
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
    mw <- tsnet::model_weights(fitobjs = object, start = start)
    weights <- mw$weights []
  } else if( n_mods > 1 & !is.null( weights ) ) {
    weights <- weights$weights[]
  } else if( n_mods == 1 ) {
    weights <- 1
  }


  ## Compute Forecasts
  # Looping over all specified models
  # I deleted the "compute_log_lik" functionality here
  forecast_obj <- lapply(fitobj, function(m){

    args <- m$arguments$fn_args
    standat <- list(Y = as.matrix(m$arguments$data),
                    "T" = m$arguments$n_t,
                    beep = args$beep,
                    prior_Beta_loc = m$arguments$priors$prior_Beta_loc,
                    prior_Beta_scale = m$arguments$priors$prior_Beta_scale,
                    prior_S = m$arguments$priors$prior_S,
                    prior_delta = m$arguments$priors$prior_delta,
                    ahead = ahead,
                    newdata = new_data,
                    compute_log_lik = compute_log_lik)

  # Obtain model name for sampling
  gqs_model <- .get_gvar_name(m)

  # backcast <- max(m$mgarchP, m$mgarchQ)
  # TODO I think this is just lag-1 here
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
  out$forecast$mean <- f_mean

  # should all samples be included?
  # if(inc_samples){
    # f_samples <- .weighted_samples(forecast_obj, weights = weights, method = method)
    # out$forecast$samples$mean <- f_samples
  # }
  # if(compute_log_lik == 1 ) {
  #   log_lik <- lapply(object.f, function(x) {
  #     rstan::extract(x, pars = "log_lik")$log_lik
  #   })
  #   out$forecast$log_lik <- log_lik
  # }


  # TODO add the arguments to the output


  class(out) <- c("tsnet_forecast", class(out))
  return(out)


}




# Helpers

# Create list of model objects
# oriented at bmgarch
tsnet_list <- function(...) {
  out <- list(...)
  class(out) <- c("tsnet_list", class(out))
  return(out)
}

# Obtain model name from stan_gvar output
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
      gqs_model <- "VAR_lkj_beep"
    }
    else{
      gqs_model <- "VAR_lkj"
    }
  }
  return(gqs_model)
}


# TODO extract weighted summaries from stan fit
# https://github.com/ph-rast/bmgarch/blob/a7d8e1c3d5a8598fa7162854802b00c25fd9a691/R/print.R#L85
# see there
