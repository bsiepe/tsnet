#' Leave-Future-Out Cross-Validation for Stan GVAR Models
#'
#' @description
#' Computes the Expected Log Predictive Density (ELPD) using Leave-One-Out (LOO) or Leave-Future-Out (LFO) approaches for models fitted using \code{\link{stan_gvar}}.
#' Currently, only the LOO method is implemented. This function is designed to work with objects of class \code{tsnet_fit}.
#'
#' @param x A \code{tsnet_fit} object obtained from \code{\link{stan_gvar}}.
#' @param type Character string specifying the type of cross-validation to perform. Options are:
#'   \itemize{
#'     \item \code{"loo"}: Leave-One-Out Cross-Validation.
#'     \item \code{"lfo"}: Leave-Future-Out Cross-Validation (not yet fully implemented).
#'   }
#'   Default is \code{"loo"}.
#' @param start Integer indicating the starting point for cross-validation. Default is 0.
#' @param ahead Integer indicating the number of steps ahead for predictions in Leave-Future-Out mode. Default is 1.
#' @param mode Character string specifying the cross-validation mode when \code{type = "lfo"}. Options are:
#'   \itemize{
#'     \item \code{"backward"}: LFO with backward refitting.
#'     \item \code{"forward"}: LFO with forward refitting (not yet implemented).
#'     \item \code{"exact"}: Exact LFO refitting at each time point (not yet implemented).
#'   }
#'   Default is \code{"backward"}.
#' @param k_thres Numeric threshold for Pareto k values. Only used in \code{lfo} mode to determine when to refit the model. Default is 0.7.
#' @param ... Additional arguments passed to underlying methods.
#'
#' @details
#' This function evaluates the predictive performance of \code{stan_gvar} models using cross-validation.
#' When \code{type = "loo"}, the function performs Leave-One-Out Cross-Validation using efficient importance sampling via the \pkg{loo} package.
#' For \code{type = "lfo"}, Leave-Future-Out Cross-Validation is conceptually supported but not yet fully implemented.
#'
#' @importFrom loo relative_eff loo psis pareto_k_values weights.importance_sampling
#'
#' @return A list with the following components:
#' \describe{
#'   \item{elpd}{Expected Log Predictive Density.}
#'   \item{mode}{Mode of cross-validation, as specified in the \code{mode} parameter.}
#' }
#'
#' @note
#' The LFO method is a work in progress and may not yet support all modes or configurations. The \code{backward} mode partially implements the LFO algorithm but does not include forecasting functionality.
#' This function borrows heavily from \code{\link[bmgarch:loo.bmgarch]{bmgarch::loo.bmgarch}}.
#'
#' @keywords internal
#' @noRd
#'

lfo.stan_gvar <- function(x,
                          type = "loo",
                          start = 0,
                          ahead = 1,
                          mode = "backward",
                          k_thres = 0.7,
                          ...) {
  obj <- x

  # Check if x is a tsnet_fit object
  if (!inherits(obj, "tsnet_fit")) {
    stop("x must be a tsnet_fit object.")
  }

  # Check if type is valid
  if (!type %in% c("loo", "lfo")) {
    stop("type must be one of 'loo' or 'lfo'.")
  }

  # Check if mode is valid
  if (!mode %in% c("backward", "forward", "exact")) {
    stop("mode must be one of 'backward', 'forward', or 'exact'.")
  }



  # Setup output list
  out <- list()
  out$type <- type
  out$mode <- mode

  # Obtain the number of observations in log likelihood (one less than full obs)
  # TODO should I change this?
  n <- obj$stan_fit@par_dims$log_lik + 1

  if (type == "loo") {

    # Obtain log-likelihood
    ll_full <- .log_lik_tsnet(obj)
    ll <- ll_full[, (start + 1):n]

    ## obtain chain id vector for relative_eff
    n_chains <- obj$stan_fit@sim$chains
    n_samples <- obj$stan_fit@sim$iter - obj$stan_fit@sim$warmup
    chain_id <- rep(seq_len(n_chains), each = n_samples)
    r_eff <- loo::relative_eff(exp(ll), chain_id = chain_id)
    backcast_loo <- loo::loo(ll, r_eff =  r_eff)
    out$elpd <- backcast_loo$estimates[, 'Estimate']['elpd_loo']
  }

  else if (type == "lfo") {
    lll <- matrix(nrow = dim(.log_lik_tsnet(obj))[1], ncol = n)
    approx_elpds_1sap <- rep(NA, n)
    ks <- NULL
    refits <- NULL

    if (mode == "backward") {
      # Start with full model, go backwards until pareto k threshold exceeded
      # then refit
      fit_past <- x
      i_refit <- n

      # loop across possible backfits
      for (i in seq((n - ahead), start, by = -ahead)) {
        # obtain log-likelihood
        ll[, (i + 1):(i + ahead)] <- .log_lik_tsnet(fit_past)[, (i + 1):(i + ahead)]
        # compute log of raw importance ratios
        # TODO check this
        logratio <- .sum_log_ratios(ll, (i + 1):i_refit)
        # compute pareto-smoothed importance sampling
        psis_obj <- suppressWarnings(loo::psis(logratio))
        # obtain k-vals
        k <- loo::pareto_k_values(psis_obj)
        ks <- c(ks, k)

        # check if k-threshold is exceeded
        if (k > k_thres) {
          # refit the model
          i_refit <- i
          refits <- c(refits, i)
          past <- 1:i
          oos <- (i + 1):(i + ahead)

          # contains data for the past
          # TODO need to check data structure (array/matrix)
          df_past <- obj$args$data[past, , drop = FALSE]
          # contains data for the out-of-sample
          # TODO why both combined?
          df_oos <- obj$args$data[c(past, oos), , drop = FALSE]
          # refitted model
          fit_past <- .refit(obj,
                             data = df_past)

          # then forecast
          # fc <- forecast(...)
          # ll[, (i+1):(i+ahead) ] <- fc$forecast$ll[[1]]
          # if(ahead == 1 ) {
          #   approx_elpds_1sap[ i+1 ] <- .log_mean_exp(ll, i+1 ])
          # } else {
          #   approx_elpds_1sap[(i+1):(i+ahead)] <-
          #     apply(ll[, (i+1):(i+ahead) ], MARGIN = 2, FUN = .log_mean_exp )
        }
        # if k threshold not exceeded
        else {
          lw <- loo:weights.importance_sampling(psis_obj, normalize = TRUE)[, 1]
          if (ahead == 1) {
            approx_elpds_1sap[i + 1] <-  .log_sum_exp(lw + ll[, i + 1])
          } else {
            approx_elpds_1sap[(i + 1):(i + ahead)] <-
              apply((lw + ll[, (i + 1):(i + ahead)]), 2, .log_sum_exp)
          }
        }
      } # end for loop
      out$elpd <- approx_elpds_1sap

    }
    else if (mode == "forward") {
      stop("Forward LFO not implemented yet.")
    }
    else if (mode == "exact") {
      stop("Exact LFO not implemented yet.")
      # refit at each time point
      k_thres <- 0
      refits <- n - (start+1)

      # todo obtain data
      df_dat <- obj$args$data
      exact_elpds_1sap <- rep(NA, n)

      # for(i in start:( n-1 ) ) {
        # fit_start <- .refit(obj,
        #                     data = df_dat[1:i, , drop = FALSE])
        # fc <- tsnet::forecast(fit_start,
        #                       ahead = ahead,
        #                       newdata = df_dat[ i+1, ,drop = FALSE])
        # loglik[, i+1 ] <-  fc$forecast$log_lik[[1]]
      # }
      # loglik_exact <- loglik[, (start+1):n]
      # exact_elpds_1sap <- apply(loglik_exact, 2, .log_mean_exp )
      # out <- exact_elpds_1sap

    }

  }


  class(out) <- c("lfo.stan_gvar", class(out))
  return(out)
}



# Function for refitting a model with given specification
# needed in lfo function above
#' DOES NOT WORK YET WITH CUSTOM PRIORS
#' @keywords internal
#' @noRd
.refit <- function(x, data) {

  args <- x$arguments$fn_args

  fit <- stan_gvar(data,
                   beep = args$beep,
                   priors = x$args$priors,
                   method = args$method,
                   cov_prior = args$cov_prior,
                   rmv_overnight = args$rmv_overnight,
                   iter_sampling = args$iter_sampling,
                   iter_warmup = args$iter_warmup,
                   n_chains = args$n_chains,
                   n_cores = args$n_cores,
                   center_only = args$center_only,
                   ahead = args$ahead,
                   compute_log_lik = args$compute_log_lik)
  return(fit)
}




# Following functions are obtained from bmgarch/lfo-cv tutorial
# Helper to obtain log-likelihood
#' @keywords internal
#' @noRd

.log_lik_tsnet <- function(x) {
  rstan::extract(x$stan_fit, pars = "log_lik")$log_lik
}

# more stable than log(sum(exp(x)))
#' @keywords internal
#' @noRd
.log_sum_exp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

## more stable than log(mean(exp(x)))
##' @keywords internal
##' @noRd
.log_mean_exp <- function(x) {
  .log_sum_exp(x) - log(length(x))
}

## compute log of raw importance ratios
## sums over observations *not* over posterior samples
##' @keywords internal
##' @noRd
.sum_log_ratios <- function(ll, ids = NULL) {
  if (!is.null(ids))
    ll <- ll[, ids , drop = FALSE]
  - rowSums(ll)
}

##' obtain relative efficiency
##' @keywords internal
##' @noRd
.rel_eff <- function(ll, x) {
  warmup <- x$stan_fit@sim$warmup
  iter <- x$stan_fit@sim$iter
  n_chains <- x$stan_fit@sim$chains
  loo::relative_eff(exp(ll), chain_id = rep(1:n_chains, each = iter - warmup))
}
