#' Leave-Future Out Cross-Validation
#'
#' @description
#' \code{lfo} computes the expected log predictive density (ELPD) with different approaches.
#' ... The function borrows heavily from \code{\link[bmgarch:loo.bmgarch]{bmgarch::loo.bmgarch}}.
#'
#'
#' @param x A \code{tsnet_fit object} obtained from \code{\link{stan_gvar}}.
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
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

  # Temporary stop because lfo not implemented yet
  if (type == "lfo") {
    stop("lfo not implemented yet.")
  }


  # Setup output list
  out <- list()
  out$type <- type
  out$mode <- mode



  if (type == "loo") {
    # Obtain the number of observations in log likelihood (one less than full obs)
    # TODO should I change this?
    n <- obj$stan_fit@par_dims$log_lik + 1

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
        ll[, (i + 1):(i + ahead)] <- .log_lik(fit_past)[, (i + 1):(i + ahead)]
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
          lw <- weights(psis_obj, normalize = TRUE)[, 1]
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
##' @keywords internal
#' DOES NOT WORK YET WITH CUSTOM PRIORS
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
                   center_only = args$center_only)
  return(fit)
}




# Following functions are obtained from bmgarch/lfo-cv tutorial

# Helper to obtain log-likelihood
.log_lik_tsnet <- function(x) {
  rstan::extract(x$stan_fit, pars = "log_lik")$log_lik
}


##' @keywords internal
## more stable than log(mean(exp(x)))
.log_mean_exp <- function(x) {
  .log_sum_exp(x) - log(length(x))
}

##' @keywords internal
## compute log of raw importance ratios
## sums over observations *not* over posterior samples
.sum_log_ratios <- function(ll, ids = NULL) {
  if (!is.null(ids))
    ll <- ll[, ids , drop = FALSE]
  - rowSums(ll)
}


##' @keywords internal
##' obtain relative efficiency
.rel_eff <- function(ll, x) {
  warmup <- x$stan_fit@sim$warmup
  iter <- x$stan_fit@sim$iter
  n_chains <- x$stan_fit@sim$chains
  loo::relative_eff(exp(ll), chain_id = rep(1:n_chains, each = iter - warmup))
}
