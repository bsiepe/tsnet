#' Prepare posterior samples for plotting
#'
#' @description This function extracts the posterior samples of partial
#' correlations and beta coefficients from a Bayesian GVAR model that was fitted
#' with [BGGM::var_estimate()] or [stan_gvar()].
#' The function is not intended to be called directly by the user,
#' but is used internally by the package.
#'
#' @param fitobj A [BGGM::var_estimate()] or [stan_gvar()]
#' fit object from which to extract the posterior samples.
#' @param burnin An integer specifying the number of initial samples to discard
#'   as burn-in. Default is 0.
#'
#' @return A matrix containing the posterior samples of partial correlations and
#'   beta coefficients. The rows represent the samples, and the columns
#'   represent the different parameters. The first set of columns represent the
#'   partial correlations, and the second set of columns represent the beta
#'   coefficients. Each column is named according to the corresponding
#'   parameter.
#'
#' @examples
#' \dontrun{
#' # Load simulated time series data
#' data(ts_data)
#' example_data <- ts_data[1:100,1:4]
#'
#' # Estimate a GVAR model
#' fit <- stan_gvar(example_data, n_chains = 2)
#'
#' # Extract posterior samples
#' plot_samples <- prepare_samples_plot (fit)
#' }
#'
#' @noRd
prepare_samples_plot <- function(fitobj,
                                   burnin = 0){

  if(inherits(fitobj, "var_estimate")) {
    # All nodes
    p <- fitobj$p
    # iterations (50 is default burnin in BGGM)
    iter <- fitobj$iter + 50

    # Number of partial corelations
    pcors_total <- p * (p - 1) * 0.5

    # identity matrix
    I_p <- diag(p)

    # pcor samples
    # 50 is default burnin in BGGM
    pcor_samples <-
      matrix(
        fitobj$fit$pcors[, , (burnin+1):(iter)][upper.tri(I_p)],
        nrow =  iter - burnin,
        ncol = pcors_total,
        byrow = TRUE
      )


    # column names
    cn <- colnames(fitobj$Y)

    if(is.null(cn)){

      col_names <- sapply(1:p, function(x)  paste(1:p, x, sep = "--"))[upper.tri(I_p)]

    } else {

      col_names <- sapply(cn, function(x)  paste(cn, x, sep = "--"))[upper.tri(I_p)]
    }
    colnames(pcor_samples) <- col_names
    posterior_samples <- pcor_samples

    n_beta_terms <- nrow(fitobj$beta_mu)
    beta_samples <- fitobj$fit$beta

    col_names <- colnames(fitobj$Y)
    beta_terms <- colnames(fitobj$X)

    beta_start <- matrix(beta_samples[1:n_beta_terms,1, (burnin+1):(iter)],
                         nrow = iter - burnin, n_beta_terms, byrow = TRUE)


    colnames(beta_start) <- paste0(col_names[1], "_",  beta_terms)

    for(i in 2:p){

      # beta next
      beta_i <- matrix(beta_samples[1:n_beta_terms, i, (burnin+1):(iter)],
                       nrow = iter - burnin,
                       n_beta_terms,
                       byrow = TRUE)

      # colnames
      colnames(beta_i) <- paste0(col_names[i], "_",  beta_terms)

      # beta combine
      beta_start <- cbind(beta_start, beta_i)

    }

    posterior_samples <-  cbind(posterior_samples, beta_start)
    return(posterior_samples)
  }

  if(inherits(fitobj, "tsnet_fit")) {
    fitobj_conv <- stan_fit_convert(fitobj,
                                    return_params = c("beta", "pcor"))

    # Number of variables
    p <- fitobj$arguments$p
    burnin <- burnin
    iter <- dim(fitobj_conv$fit$beta)[3]

    # Number of partial correlations
    pcors_total <- p * (p-1) * 0.5

    # Identity
    I_p <- diag(p)

    # Get the samples
    pcor_samples <- matrix(fitobj_conv$fit$pcors[,,][upper.tri(I_p)],
                           nrow = dim(fitobj_conv$fit$pcors)[3],
                           ncol = pcors_total,
                           byrow = TRUE)

    # column names
    cnames <- fitobj$arguments$cnames

    if(is.null(cnames)){

      col_names <- sapply(1:p, function(x)  paste(1:p, x, sep = "--"))[upper.tri(I_p)]

    } else {

      col_names <- sapply(cnames, function(x)  paste(cnames, x, sep = "--"))[upper.tri(I_p)]
    }

    colnames(pcor_samples) <- col_names
    posterior_samples <- pcor_samples

    n_beta_terms <- nrow(fitobj_conv$beta_mu)
    beta_samples <- fitobj_conv$fit$beta

    col_names <- colnames(fitobj_conv$beta_mu)
    beta_terms <- rownames(fitobj_conv$beta_mu)

    beta_start <- matrix(beta_samples[1:n_beta_terms,1, (burnin+1):(iter)],
                         nrow = iter - burnin, n_beta_terms, byrow = TRUE)


    colnames(beta_start) <- paste0(col_names[1], "_",  beta_terms)

    for(i in 2:p){

      # beta next
      beta_i <- matrix(beta_samples[1:n_beta_terms, i, (burnin+1):(iter)],
                       nrow = iter - burnin,
                       n_beta_terms,
                       byrow = TRUE)

      # colnames
      colnames(beta_i) <- paste0(col_names[i], "_",  beta_terms)

      # beta combine
      beta_start <- cbind(beta_start, beta_i)

    }

    samps <-  cbind(posterior_samples, beta_start)
  }
  return(samps)
}


#' Convert Stan Fit to Array of Samples
#'
#' This function converts a Stan fit object into an array of samples for the
#' temporal coefficients and the innovation covariance or partial correlation
#' matrices. It supports rstan as a backend. It can be used to convert models
#' fit using [stan_gvar()] into 3D arrays, which is the standard data structure
#' used in `tsnet`. The function allows to select which parameters should be
#' returned.
#'
#' @param stan_fit A Stan fit object obtained from rstan or a tsnet_fit object
#'   from [stan_gvar()].
#' @param return_params A character vector specifying which parameters to
#'   return. Options are "beta" (temporal network), "sigma" (innovation
#'   covariance), and "pcor" (partial correlations). Default is
#'   c("beta","sigma", "pcor").
#'
#' @return A list containing 3D arrays for the selected parameters. Each array
#'   represents the posterior samples for a parameter, and each slice of the
#'   array represents a single iteration.
#' @importFrom rstan extract
#' @importFrom posterior as_draws_matrix
#' @examples
#' \dontrun{
#' data(ts_data)
#' example_data <- ts_data[1:100,1:3]
#' fit <- stan_gvar(data = example_data,
#'                  n_chains = 2,
#'                  n_cores = 1)
#' samples <- stan_fit_convert(fit, return_params = c("beta", "pcor"))
#' }
#'
#' @export
stan_fit_convert <- function(stan_fit,
                             return_params = c("beta", "sigma", "pcor")) {


  if (inherits(stan_fit, "tsnet_fit")) {
    stan_obj <- stan_fit$stan_fit
    cnames <- stan_fit$arguments$cnames
  } else {
    stan_obj <- stan_fit
    cnames <- NULL
  }

  # Number of variables
  p <- stan_fit$stan_fit@par_dims$Beta[1]

  if(is.null(cnames)){
    cnames <- paste0("V", 1:p)
  }
  names0 <- cnames
  names1 <- paste0(cnames, ".l1")

  # check fitting backend
  c <- class(stan_obj)

  if (attr(c, "package") == "rstan") {
    if ("beta" %in% return_params) {
      draws_beta <-
        posterior::as_draws_matrix(rstan::extract(stan_obj, pars = "Beta", permuted = FALSE))
    }
    if ("sigma" %in% return_params) {
      draws_sigma <-
        posterior::as_draws_matrix(rstan::extract(stan_obj, pars = "Sigma", permuted = FALSE))
    }
    if ("pcor" %in% return_params) {
      draws_pcor <-
        posterior::as_draws_matrix(rstan::extract(stan_obj, pars = "Rho", permuted = FALSE))
    }
  } else {

    stop("Only rstan backend within `stan_gvar` is supported at the moment.")
  }

  # Convert to array of p x p matrices
  nvar <- sqrt(ncol(draws_beta))

  return_list <- list()
  return_list$fit <- list()

  if ("beta" %in% return_params) {
    split_beta <- split(draws_beta, seq(nrow(draws_beta)))
    beta_l <- lapply(split_beta, function(x) {
      matrix(x,
             nrow = nvar,
             ncol = nvar,
             byrow = TRUE)
    })
    return_list$fit$beta <-
      array(unlist(beta_l), dim = c(nvar, nvar, nrow(draws_beta)))
    return_list$beta_mu <- apply(return_list$fit$beta, c(1, 2), mean)
    rownames(return_list$beta_mu) <- names1
    colnames(return_list$beta_mu) <- names0
  }

  if ("sigma" %in% return_params) {
    split_sigma <- split(draws_sigma, seq(nrow(draws_sigma)))
    sigma_l <- lapply(split_sigma, function(x) {
      matrix(x,
             nrow = nvar,
             ncol = nvar,
             byrow = TRUE)
    })
    return_list$fit$sigma <-
      array(unlist(sigma_l), dim = c(nvar, nvar, nrow(draws_sigma)))
    return_list$sigma_mu <- apply(return_list$fit$sigma, c(1, 2), mean)
    rownames(return_list$sigma_mu) <- cnames
    colnames(return_list$sigma_mu) <- cnames
  }

  if ("pcor" %in% return_params) {
    split_pcor <- split(draws_pcor, seq(nrow(draws_pcor)))
    pcor_l <- lapply(split_pcor, function(x) {
      matrix(x,
             nrow = nvar,
             ncol = nvar,
             byrow = TRUE)
    })
    return_list$fit$pcors <-
      array(unlist(pcor_l), dim = c(nvar, nvar, nrow(draws_pcor)))
    return_list$pcor_mu <- apply(return_list$fit$pcors, c(1, 2), mean)
    rownames(return_list$pcor_mu) <- cnames
    colnames(return_list$pcor_mu) <- cnames
  }


  class(return_list) <- c("tsnet_samples", class(return_list))

  return(return_list)

}

