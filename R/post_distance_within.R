#' Calculates distances between pairs of posterior samples using the posterior
#' samples or posterior predictive draws
#'
#' @description This function computes distances between posterior samples of a
#'   single fitted GVAR model. Thereby, it calculates the uncertainty contained
#'   in the posterior distribution, which can be used as a reference to compare
#'   two modes. Distances can be obtained either from posterior samples or
#'   posterior predictive draws. The distance between two models can currently
#'   be calculated based on three options: Frobenius norm, maximum difference,
#'   or L1 norm. Used within [compare_gvar()]. The function is not intended to
#'   be used directly by the user.
#'
#' @inheritParams compare_gvar
#' @param fitobj Fitted model object. This can be a tsnet_fit object (obtained
#'   from [stan_gvar()]), a BGGM object (obtained from [BGGM::var_estimate()]),
#'   or extracted posterior samples (obtained from [stan_fit_convert()).
#' @param pred A logical indicating whether the input is posterior predictive
#'   draws (TRUE) or posterior samples (FALSE). Default: FALSE
#' @return A list of distances between the specified pairs of fitted models. The
#'   list has length equal to the specified number of random pairs. Each list
#'   element contains two distance values, one for beta coefficients and one for
#'   partial correlations.
#' @examples
#' data(fit_data)
#' post_distance_within(fitobj = fit_data[[1]],
#' comp = "frob",
#' pred = FALSE,
#' n_draws = 100)
#'
#' @export

post_distance_within <- function(fitobj,
                                 comp,
                                 pred,
                                 n_draws = 1000,
                                 sampling_method = "random",
                                 indices = NULL,
                                 burnin = 0) {
  # storage
  dist_out <- list()


  # Helper function to only use upper triangle of matrix
  ut <- function(x) {
    matrix(x[upper.tri(x, diag = FALSE)])
  }

  # Helper to compute metric
  compute_metric <- function(a, b, metric, indices = NULL) {
    tryCatch(
      {
        if (!is.null(indices)) {
          a <- a[indices]
          b <- b[indices]
        }
        if (metric == "frob") {
          norm(a - b, type = "F")
        } else if (metric == "maxdiff") {
          max(abs(a - b))
        } else if (metric == "l1") {
          sum(abs(a - b))
        }
      },
      error = function(e) NA
    )
  }

  # define the distance function based on comp
  # draw from all posterior samples
  # currently uses identical function calls and applies different arguments below
  distance_fn_beta <- function(a, b, comp, indices) compute_metric(a, b, comp, indices)
  distance_fn_pcor <- function(a, b, comp, indices) compute_metric(a, b, comp, indices)

  # for posterior predictive approach
  if (isTRUE(pred)) {
    # Obtain number of models
    n_mod <- length(fitobj$fit)
  }

  # for posteriors of empirical models
  if (isFALSE(pred)) {
    if (!is.list(fitobj)) {
      beta_distance <- NA
      pcor_distance <- NA
      stop("Not a list.")
    } else{
    # Obtain number of posterior samples
    n_mod <- dim(fitobj$fit$beta)[3]
    }
  }

  ## Draw two random models

  # "samples" would be more fitting here than "models"
  # "model" is still a residue from posterior predictive approach
  # Draw models spaced apart so that we don't have autocorrelation from sampling
  mod_pairs <- array(NA, dim = c(2, n_draws))


  if(sampling_method == "sequential"){
    # draw from first half of samples
    mod_pairs[1, 1:(n_draws)] <- seq(51, n_draws + 50, by = 1)

    # draw from second half of samples
    mod_pairs[2, 1:(n_draws)] <- seq((n_mod / 2) + 51,
                                     (n_mod / 2) + 50 + (n_draws), by = 1)

  }
  if(sampling_method == "random"){
    # Determine the valid range of values for drawing pairs
    # Leave out middle 100 to keep certain distance between halves
    # take out burnin-in samples
    valid_range <- setdiff((burnin+1):n_mod,
                           ((n_mod/2 - 50 + (burnin/2)):(n_mod/2 + 49 + burnin/2)))

    # Draw pairs randomly from the first half and second half
    mod_pairs[1, 1:n_draws] <- sample(valid_range[1:(length(valid_range)/2)],
                                      n_draws)
    mod_pairs[2, 1:n_draws] <- sample(valid_range[(length(valid_range)/2 + 1):length(valid_range)], n_draws)


  }

  for (i in seq(n_draws)) {
    # storage
    dist_out[[i]] <- list()
    mod_one <- mod_pairs[1, i]
    mod_two <- mod_pairs[2, i]

    # if mod_one and mod_two are equal, redraw
    if (mod_one == mod_two) {
      mod_two <- sample(1:n_mod, size = 1)
    }


  if (isTRUE(pred)) {
    # conditional: ut() is only called if indices is NULL to keep correct indexing for user-specified indices
      beta_distance <- if(is.null(indices)) {
        distance_fn_beta(fitobj$fit[[mod_one]]$beta_mu,
                         fitobj$fit[[mod_two]]$beta_mu,
                         comp,
                         indices)
      } else {
        distance_fn_beta(fitobj$fit[[mod_one]]$beta_mu,
                         fitobj$fit[[mod_two]]$beta_mu,
                         comp,
                          indices$beta)
        }

      pcor_distance <- if (is.null(indices)) {
        distance_fn_pcor(ut(fitobj$fit[[mod_one]]$pcor_mu),
                         ut(fitobj$fit[[mod_two]]$pcor_mu),
                         comp,
                         indices)
      } else {
        distance_fn_pcor(fitobj$fit[[mod_one]]$pcor_mu,
                         fitobj$fit[[mod_two]]$pcor_mu,
                         comp,
                         indices$pcor)
      }

  }

  if (isFALSE(pred)) {
    # if both elements are lists
    # only use indices if probided

      beta_distance <- if (is.null(indices)){
        distance_fn_beta(fitobj$fit$beta[, , mod_one],
                         fitobj$fit$beta[, , mod_two],
                         comp,
                         indices)
      } else {
        distance_fn_beta(fitobj$fit$beta[, , mod_one],
                         fitobj$fit$beta[, , mod_two],
                         comp,
                         indices$beta)
      }
      pcor_distance <- if (is.null(indices)) {
        distance_fn_pcor(ut(fitobj$fit$pcors[, , mod_one]),
                         ut(fitobj$fit$pcors[, , mod_two]),
                         comp,
                         indices)
      } else {
        distance_fn_pcor(fitobj$fit$pcors[, , mod_one],
                         fitobj$fit$pcors[, , mod_two],
                         comp,
                         indices$pcor)
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
