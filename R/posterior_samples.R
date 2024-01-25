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
#' fit <- var_estimate (...) # TODO replace with actual fitting function
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





