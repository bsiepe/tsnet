#' Extract Posterior Samples from a BGGM Fit Object
#'
#' This function extracts the posterior samples of partial correlations and beta coefficients from a fitted BGGM model.
#'
#' @param fitobj A BGGM fit object from which to extract the posterior samples.
#' @param burnin An integer specifying the number of initial samples to discard as burn-in. Default is 500.
#'
#' @return A matrix containing the posterior samples of partial correlations and beta coefficients. The rows represent the samples, and the columns represent the different parameters. The first set of columns represent the partial correlations, and the second set of columns represent the beta coefficients. Each column is named according to the corresponding parameter.
#'
#' @examples
#' \dontrun{
#' fit <- var_estimate (...) # replace with actual fitting function
#' posterior_samples <- posterior_samples_bggm(fit)
#' }
#'
#' @export
posterior_samples_bggm <- function(fitobj,
                                   burnin = 500){

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



#' Convert Draws Matrix to List of Matrices
#'
#' This helper function transforms a matrix of draws into a list of matrices,
#' each representing a single iteration.
#'
#' @param draws_matrix A matrix of draws where each row represents an iteration.
#' @return A list of matrices, each representing a single iteration.
draws_matrix2list <- function(draws_matrix) {
  iterations_list <-
    lapply(
      X = 1:nrow(draws_matrix),
      FUN = function(X) {
        matrix(draws_matrix[X,], ncol = sqrt(ncol(draws_matrix)), byrow = FALSE)
      }
    )
  return(iterations_list)
}

#' Convert Draws Matrix to Array
#'
#' This helper function transforms a matrix of draws into a 3D array.
#'
#' @param draws_matrix A matrix of draws where each row represents an iteration.
#' @return A 3D array where each slice represents an iteration.
draws_matrix2array <- function(draws_matrix) {
  array <-
    array(t(draws_matrix),
          dim = c(sqrt(ncol(draws_matrix)),
                  sqrt(ncol(draws_matrix)),
                  nrow(draws_matrix)))

  return(array)
}

#' Convert Draws Array to Matrix
#'
#' This helper function transforms a 3D array of draws into a matrix.
#' It also allows for the removal of warmup samples.
#'
#' @param array_3d A 3D array of draws where each slice represents an iteration.
#' @param warmup An integer specifying the number of initial samples to discard as warm-up. Default is 0.
#' @return A matrix where each row represents an iteration.
draws_array2matrix <- function(array_3d,
                               warmup = 0) { # set to zero to keep everything
  iterations_list <-
    lapply(
      X = (warmup+1):dim(array_3d)[3],
      FUN = function(X) {
        as.vector(array_3d[, , X])
      }
    )
  matrix <- do.call(rbind, iterations_list)
  return(matrix)
}

#' Convert Stan Fit to Array of Samples
#'
#' This function converts a Stan fit object into an array of samples for each parameter (Beta, Sigma, Pcor).
#' It supports rstan as a backend. The function allows to select which parameters should be returned.
#'
#' @param stan_fit A Stan fit object obtained from rstan.
#' @param return_params A character vector specifying which parameters to return. Options are "beta", "sigma", and "pcor". Default is c("beta", "sigma", "pcor").
#'
#' @return A list containing 3D arrays for the selected parameters. Each array represents the posterior samples for a parameter, and each slice of the array represents a single iteration.
#'
#' @examples
#' \dontrun{
#' samples <- stan_fit_convert(stan_fit, return_params = c("beta", "pcor"))
#' }
#'
#' @export
stan_fit_convert <- function(stan_fit, 
                             return_params = c("beta", "sigma", "pcor")) {
  # check fitting backend
  c <- class(stan_fit)

  if (attr(c, "package") == "rstan") {
    if ("beta" %in% return_params) {
      draws_beta <- posterior::as_draws_matrix(rstan::extract(stan_fit, pars = "Beta", permuted = FALSE))
    }
    if ("sigma" %in% return_params) {
      draws_sigma <- posterior::as_draws_matrix(rstan::extract(stan_fit, pars = "Sigma", permuted = FALSE))
    }
    if ("pcor" %in% return_params) {
      draws_pcor <- posterior::as_draws_matrix(rstan::extract(stan_fit, pars = "Rho", permuted = FALSE))
    }
  } else {
    stop("Only rstan backend within `stan_gvar` is supported at the moment.")
  }

  # Convert to array of p x p matrices
  nvar <- sqrt(ncol(draws_beta))

  return_list <- list()

  if ("beta" %in% return_params) {
    split_beta <- split(draws_beta, seq(nrow(draws_beta)))
    beta_l <- lapply(split_beta, function(x) {
      matrix(x, 
            nrow = nvar, 
            ncol = nvar, 
            byrow = TRUE)
    })
    return_list$beta <- array(unlist(beta_l), dim = c(nvar, nvar, nrow(draws_beta)))
  }

  if ("sigma" %in% return_params) {
    split_sigma <- split(draws_sigma, seq(nrow(draws_sigma)))
    sigma_l <- lapply(split_sigma, function(x) {
      matrix(x, 
            nrow = nvar, 
            ncol = nvar, 
            byrow = TRUE)
    })
    return_list$sigma <- array(unlist(sigma_l), dim = c(nvar, nvar, nrow(draws_sigma)))
  }

  if ("pcor" %in% return_params) {
    split_pcor <- split(draws_pcor, seq(nrow(draws_pcor)))
    pcor_l <- lapply(split_pcor, function(x) {
      matrix(x, 
            nrow = nvar, 
            ncol = nvar, 
            byrow = TRUE)
    })
    return_list$pcor <- array(unlist(pcor_l), dim = c(nvar, nvar, nrow(draws_pcor)))
  }

  return(return_list)
}

