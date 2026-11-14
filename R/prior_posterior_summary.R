#' Compare Prior and Posterior Distributions of a GVAR Model
#'
#' @description Summarizes, for every edge of the temporal and the
#'   contemporaneous network, the mean and standard deviation of the prior
#'   and the posterior distribution, as well as the posterior contraction
#'   \eqn{1 - Var_{post} / Var_{prior}}. A contraction close to 1 indicates
#'   that the data were much more informative than the prior, whereas a
#'   contraction close to 0 indicates that the posterior is dominated by the
#'   prior. Negative values occur if the posterior is wider than the prior.
#'
#' @param fit A \code{tsnet_fit} object obtained from \code{\link{stan_gvar}}
#'   with \code{return_all = TRUE} (the default).
#' @param prior_fit Optional \code{tsnet_fit} object obtained from
#'   \code{\link{stan_gvar}} with \code{prior_only = TRUE}, using the same
#'   priors and \code{cov_prior} as \code{fit}. Required to summarize the
#'   prior of the partial correlations, which is not available in closed
#'   form in general. Default: \code{NULL}.
#'
#' @details Prior means and standard deviations of the temporal coefficients
#'   are taken from the normal priors stored in \code{fit}. Those of the
#'   partial correlations are estimated from the draws of \code{prior_fit}
#'   and are \code{NA} if no prior-only fit is provided.
#'
#' @return A data frame with one row per edge and the columns \code{type}
#'   (\code{"temporal"} or \code{"contemporaneous"}), \code{from}, \code{to},
#'   \code{prior_mean}, \code{prior_sd}, \code{post_mean}, \code{post_sd},
#'   and \code{contraction}. For temporal edges, \code{from} is the lagged
#'   predictor and \code{to} is the outcome. Contemporaneous edges are
#'   undirected and listed once.
#'
#' @seealso \code{\link{stan_gvar}}, \code{\link{gvar_priors}}
#'
#' @examples
#' \donttest{
#' data(ts_data)
#' example_data <- ts_data[1:100, 1:3]
#' fit <- stan_gvar(example_data,
#'   cov_prior = "IW",
#'   n_chains = 2
#' )
#' prior_fit <- stan_gvar(example_data,
#'   cov_prior = "IW",
#'   n_chains = 2,
#'   prior_only = TRUE
#' )
#' prior_posterior_summary(fit, prior_fit)
#' }
#'
#' @export
prior_posterior_summary <- function(fit, prior_fit = NULL) {
  if (!inherits(fit, "tsnet_fit") || isTRUE(fit$arguments$prior_only)) {
    stop("'fit' must be a 'tsnet_fit' object that is not a prior-only fit.")
  }
  priors <- fit$arguments$priors
  if (is.null(priors)) {
    stop("'fit' does not contain its priors. Refit with 'return_all = TRUE'.")
  }
  K <- fit$arguments$p
  cnames <- fit$arguments$cnames
  if (is.null(cnames)) {
    cnames <- paste0("V", seq_len(K))
  }

  beta_draws <- .extract_draws(fit, "Beta")
  # rows of Beta are outcomes, columns are lagged predictors
  temporal <- data.frame(
    type = "temporal",
    from = cnames[as.vector(col(diag(K)))],
    to = cnames[as.vector(row(diag(K)))],
    prior_mean = as.vector(priors$prior_Beta_loc),
    prior_sd = as.vector(priors$prior_Beta_scale),
    post_mean = as.vector(apply(beta_draws, c(2, 3), mean)),
    post_sd = as.vector(apply(beta_draws, c(2, 3), stats::sd))
  )

  upper <- upper.tri(diag(K))
  rho_draws <- .extract_draws(fit, "Rho")
  # no contemporaneous edges for a single variable
  n_pcor <- sum(upper)
  contemporaneous <- data.frame(
    type = rep("contemporaneous", n_pcor),
    from = cnames[row(diag(K))[upper]],
    to = cnames[col(diag(K))[upper]],
    prior_mean = rep(NA_real_, n_pcor),
    prior_sd = rep(NA_real_, n_pcor),
    post_mean = apply(rho_draws, c(2, 3), mean)[upper],
    post_sd = apply(rho_draws, c(2, 3), stats::sd)[upper]
  )

  if (!is.null(prior_fit)) {
    if (!inherits(prior_fit, "tsnet_fit") ||
      !isTRUE(prior_fit$arguments$prior_only)) {
      stop("'prior_fit' must be obtained from stan_gvar(prior_only = TRUE).")
    }
    same_prior <- identical(prior_fit$arguments$priors, priors) &&
      identical(.cov_model(prior_fit), .cov_model(fit))
    if (!same_prior) {
      warning("'fit' and 'prior_fit' do not seem to use the same priors.")
    }
    prior_rho <- .extract_draws(prior_fit, "Rho")
    if (dim(prior_rho)[2] != K) {
      stop("'fit' and 'prior_fit' differ in the number of variables.")
    }
    contemporaneous$prior_mean <-
      apply(prior_rho, c(2, 3), mean, na.rm = TRUE)[upper]
    contemporaneous$prior_sd <-
      apply(prior_rho, c(2, 3), stats::sd, na.rm = TRUE)[upper]
  }

  out <- rbind(temporal, contemporaneous)
  out$contraction <- 1 - out$post_sd^2 / out$prior_sd^2
  out
}

# stan model without the suffix of the models that remove overnight effects
.cov_model <- function(fit) {
  sub("_beep$", "", fit$stan_fit@model_name)
}
