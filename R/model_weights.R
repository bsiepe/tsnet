#' Compute Model Weights Based on Leave-Future-Out Cross-Validation
#'
#' @description
#' This function computes model weights based on Leave-Future-Out (LFO) cross-validation
#' for a list of fitted GVAR models. It leverages the LFO method to evaluate predictive performance
#' and calculates weights using the \code{\link[loo:loo_model_weights]{loo::loo_model_weights}} function.
#'
#' @param fits A list of \code{tsnet_fit} objects or an object of class \code{tsnet_list}.
#' @inheritParams lfo.stan_gvar
#' @param return_lfo Logical. If \code{TRUE}, the function also returns the LFO objects. Default is \code{FALSE}.
#' @param ... Additional arguments (currently not used).
#'
#' @details
#' The function processes the input list of fitted models, ensuring compatibility with the \code{tsnet_list} class.
#' It computes the LFO cross-validation results for each model in the list using \code{\link{lfo.stan_gvar}}
#' and extracts log-likelihoods and relative efficiency information. The final model weights are then calculated
#' based on the log-likelihoods and relative efficiency values using the \code{\link[loo:loo_model_weights]{loo::loo_model_weights}} function.
#'
#' The LFO cross-validation options (\code{start}, \code{ahead}, \code{mode}, \code{k_thres})
#' are inherited from \code{\link{lfo.stan_gvar}}.
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
#'
#' @examples
#' \dontrun{
#' # Assuming `fit1` and `fit2` are tsnet_fit objects
#' fits <- tsnet_list(fit1, fit2)
#' weights <- model_weights(fits)
#' print(weights)
#' }
#' @keywords internal
#' @noRd

model_weights <- function(fits,
                          start = 0,
                          ahead = 1,
                          mode = "backward",
                          k_thres = 0.7,
                          return_lfo = FALSE,
                          ...) {

  # browser()
  # Check if input is a list
  if (!is.list(fits)) {
    stop("Input must be a list of tsnet_fit objects")
  }

  # Check if input is tsnet_list
  # if not, convert to tsnet_list
  if (!inherits(fits, "tsnet_list")) {
    fits <- tsnet_list(fits)
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
  loo_method = "stacking"

  weights <- loo::loo_model_weights(
    ll_list,
    method = loo_method,
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
