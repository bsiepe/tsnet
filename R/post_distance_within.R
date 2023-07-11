#' Calculates distances between pairs of posterior samples using the posterior samples or posterior predictive draws
#'
#' This function computes distances between posterior samples of a fitted gVAR model.
#' Distances can be obtained either from posterior samples or posterior predictive draws.
#' The distance between two models can currently be calculated based on three options: Frobenius norm, maximum difference, or L1 norm.
#' Used within [compare_gvar()].
#'
#' @param fitobj A BGGM var_estimate fit object.
#' @param comp A character string indicating the type of distance between models that should be calculated. The options include: "frob" (Frobenius norm), "maxdiff" (maximum difference), or "l1" (L1 norm).
#' @param pred A logical indicating whether the input is posterior predictive draws (TRUE) or posterior samples (FALSE).
#' @param draws An integer specifying the number of random pairs of models that should be compared.
#'
#' @return A list of distances between the specified pairs of fitted models. The list has length equal to the specified number of random pairs. Each list element contains two distance values, one for beta coefficients and one for partial correlations.
#'
#'
#' @export post_distance_within
#'

post_distance_within <- function(fitobj,
                                 comp,
                                 pred, # posterior predictive?
                                 draws = 1000) {
  # storage
  dist_out <- list()


  # for posterior predictive approach
  if (isTRUE(pred)) {
    # define the distance function based on comp
    distance_fn_beta <- switch(comp,
      frob = {
        function(x, y, mod_one, mod_two) norm(x$fit[[mod_one]]$beta_mu - y$fit[[mod_two]]$beta_mu, type = "F")
      },
      maxdiff = {
        function(x, y, mod_one, mod_two) max(abs((x$fit[[mod_one]]$beta_mu - y$fit[[mod_two]]$beta_mu)))
      },
      l1 = {
        function(x, y, mod_one, mod_two) sum(abs((x$fit[[mod_one]]$beta_mu - y$fit[[mod_two]]$beta_mu)))
      }
    )
    distance_fn_pcor <- switch(comp,
      frob = {
        function(x, y, mod_one, mod_two) norm(x$fit[[mod_one]]$pcor_mu - y$fit[[mod_two]]$pcor_mu, type = "F")
      },
      maxdiff = {
        function(x, y, mod_one, mod_two) max(abs((x$fit[[mod_one]]$pcor_mu - y$fit[[mod_two]]$pcor_mu)))
      },
      l1 = {
        function(x, y, mod_one, mod_two) sum(abs((x$fit[[mod_one]]$pcor_mu - y$fit[[mod_two]]$pcor_mu)))
      }
    )

    # Obtain number of models
    n_mod <- length(fitobj$fit)
  }


  # for posteriors of empirical models
  if (isFALSE(pred)) {
    # define the distance function based on comp
    # draw from all posterior samples

    distance_fn_beta <- switch(comp,
      frob = {
        function(x, y, mod_one, mod_two) norm(x$fit$beta[, , mod_one] - y$fit$beta[, , mod_two], type = "F")
      },
      maxdiff = {
        function(x, y, mod_one, mod_two) max(abs((x$fit$beta[, , mod_one] - y$fit$beta[, , mod_two])))
      },
      l1 = {
        function(x, y, mod_one, mod_two) sum(abs((x$fit$beta[, , mod_one] - y$fit$beta[, , mod_two])))
      }
    )
    distance_fn_pcor <- switch(comp,
      frob = {
        function(x, y, mod_one, mod_two) norm(x$fit$pcors[, , mod_one] - y$fit$pcors[, , mod_two], type = "F")
      },
      maxdiff = {
        function(x, y, mod_one, mod_two) max(abs((x$fit$pcors[, , mod_one] - y$fit$pcors[, , mod_two])))
      },
      l1 = {
        function(x, y, mod_one, mod_two) sum(abs((x$fit$pcors[, , mod_one] - y$fit$pcors[, , mod_two])))
      }
    )

    # Obtain number of posterior samples
    n_mod <- dim(fitobj$fit$beta)[3]
  }


  ## Draw two random models
  # delete burn-in iterations (to 50)
  # n_mod <- n_mod[-c(1:50)]

  # "samples" would be more fitting here than "models"
  # "model" is still a residue from posterior predictive approach
  # Draw models spaced apart so that we don't have autocorrelation from sampling
  mod_pairs <- array(NA, dim = c(2, draws))
  # draw from first half of samples
  mod_pairs[1, 1:(draws)] <- seq(51, draws + 50, by = 1)

  # draw from second half of samples
  mod_pairs[2, 1:(draws)] <- seq((n_mod / 2) + 51, (n_mod / 2) + 50 + (draws), by = 1)

  # mod_pairs <- replicate(draws, sample(1:n_mod, size = 2, replace = TRUE))

  for (i in seq(draws)) {
    # storage
    dist_out[[i]] <- list()
    mod_one <- mod_pairs[1, i]
    mod_two <- mod_pairs[2, i]

    # if mod_one and mod_two are equal, redraw
    if (mod_one == mod_two) {
      mod_two <- sample(1:n_mod, size = 1)
    }

    ## Check if estimation worked
    # Should be unneccessary if non-converged attempts were deleted
    if (isTRUE(pred)) {
      if (!is.list(fitobj$fit[[mod_one]]) | !is.list(fitobj$fit[[mod_two]])) {
        beta_distance <- NA
        pcor_distance <- NA
        stop("Not a list.")
      }
      # if both elements are lists
      else {
        beta_distance <- distance_fn_beta(fitobj, fitobj, mod_one, mod_two)
        pcor_distance <- distance_fn_pcor(fitobj, fitobj, mod_one, mod_two)
      }
    }

    if (isFALSE(pred)) {
      if (!is.list(fitobj) | !is.list(fitobj)) {
        beta_distance <- NA
        pcor_distance <- NA

        stop("Not a list.")
      }
      # if both elements are lists
      else {
        beta_distance <- distance_fn_beta(fitobj, fitobj, mod_one, mod_two)
        pcor_distance <- distance_fn_pcor(fitobj, fitobj, mod_one, mod_two)
      }
    }

    # Store results
    dist_out[[i]]$comp <- comp
    dist_out[[i]]$mod_one <- mod_one
    dist_out[[i]]$mod_two <- mod_two
    dist_out[[i]]$beta <- beta_distance
    dist_out[[i]]$pcor <- pcor_distance
  } # end for loop
  out <- do.call(rbind.data.frame, dist_out)


  return(out)
}
