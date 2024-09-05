#' Compute model weights based on LFO
#'
#' This function computes model weights based on Leave-Future-Out (LFO) cross-validation.
#' It takes a list of fitted GVAR models and calculates the model weights using the LFO method.
#' This function is oriented at the \code{\link[bmgarch:model_weights]{bmgarch::model_weights}} function.
#'
#' @param fits A list of \code{tsnet_fit} objects or an object of class \code{tsnet_list}.
#' @inheritParams lfo.stan_gvar
#' @param return_lfo Logical. If TRUE, the function also returns the LFO objects. Default is FALSE.
#' @param ... Currently not in use
#'
#' @details The function first checks if the input is a list of \code{tsnet_fit} objects.
#' If not, it stops with an error. If the input is not of class \code{tsnet_list}, it converts it to \code{tsnet_list}.
#' It then computes the LFO for each model in the list and extracts the log-likelihoods.
#' The relative efficiency information is also obtained for each model.
#' Finally, the model weights are computed using the \code{\link[loo:loo_model_weights]{loo::model_weights}} function.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{weights}{The computed model weights.}
#'   \item{ll_list}{A list of log-likelihoods for each model.}
#'   \item{r_eff_list}{A list of relative efficiency information for each model.}
#'   \item{lfo_list}{(Optional) A list of LFO objects if \code{return_lfo} is \code{TRUE}.}
#' }
#'
#' @importFrom loo loo_model_weights
#' @importFrom tsnet tsnet_list
#'
#' @examples
#' \donttest{
#' # Assuming `fit1` and `fit2` are tsnet_fit objects
#' fits <- tsnet_list(fit1, fit2)
#' weights <- model_weights(fits)
#' print(weights)
#' }
#' @export
model_weights <- function(fits,
                          start = 0,
                          ahead = 1,
                          mode = "backward",
                          k_thres = 0.7,
                          return_lfo = FALSE,
                          ...) {

  browser()
  # Check if input is a list
  if (!is.list(fits)) {
    stop("Input must be a list of tsnet_fit objects")
  }

  # Check if input is tsnet_list
  # if not, convert to tsnet_list
  if (!inherits(fits, "tsnet_list")) {
    fits <- tsnet::tsnet_list(fits)
  }

  ## Compute model weights based on LFO
  lfo_list <- lapply(
    fits,
    FUN = lfo.stan_gvar,
    type = "lfo",
    start = start,
    ahead = ahead,
    mode = mode
  )

  # TODO this is not correct yet
  # we need to extract the log-likelihoods from the LFO objects, not the elpd
  # see https://github.com/ph-rast/bmgarch/blob/master/R/model_weights.R#L3
  ll_list <- lapply(
    lfo_list,
    FUN = function(x)
      x$loglik
  )

  # Obtain rel_eff info
  r_eff_list <- list()
  for (i in seq_along(lfo_list)) {
    r_eff_list[[i]] <- .rel_eff(ll_list[[i]], fits[[i]])
  }


  # Compute model weights
  # adapted from bmgarch::model_weights
  weights <- loo::loo_model_weights(
    ll_list,
    method = method,
    r_eff_list = r_eff_list,
    optim_control = list(reltol = 1e-10)
  )

  out <- list()
  out$weights <- weights
  out$ll_list <- ll_list
  out$r_eff_list <- r_eff_list

  if (isTRUE(return_lfo)) {
    out$lfo_list <- lfo_list
  }
  class(out) <- c("model_weights", class(out))
  return(out)

}
