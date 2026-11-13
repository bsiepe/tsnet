#' Specify Priors for a Bayesian GVAR Model
#'
#' @description Returns the default priors of \code{\link{stan_gvar}} for a
#'   model with \code{K} variables, optionally with some of them replaced.
#'   In contrast to passing an incomplete list to the \code{priors} argument of
#'   \code{\link{stan_gvar}}, misspelled names and priors with wrong
#'   dimensions result in an error.
#'
#' @param K Number of variables in the model.
#' @param ... Named priors that replace the defaults. Possible names are
#'   \code{prior_Beta_loc}, \code{prior_Beta_scale}, \code{prior_Rho_loc},
#'   \code{prior_Rho_scale}, \code{prior_Rho_marginal}, \code{prior_S}, and
#'   \code{prior_Eta}. See the section on priors in \code{\link{stan_gvar}}
#'   for their meaning. For matrices of temporal coefficients, rows are
#'   outcomes and columns are lagged predictors.
#'
#' @return A named list of priors that can be passed to the \code{priors}
#'   argument of \code{\link{stan_gvar}}.
#'
#' @seealso \code{\link{stan_gvar}}, \code{\link{priors_from_fit}}
#'
#' @examples
#' # default priors for three variables
#' str(gvar_priors(K = 3))
#'
#' # informative prior on the autoregressive effects
#' gvar_priors(K = 3, prior_Beta_loc = diag(0.3, 3))$prior_Beta_loc
#'
#' @export
gvar_priors <- function(K, ...) {
  if (!.is_count(K) || K < 1) {
    stop("'K' must be a single positive integer.")
  }
  defaults <- .default_priors(K)
  overrides <- list(...)
  if (length(overrides) > 0 &&
    (is.null(names(overrides)) || any(names(overrides) == ""))) {
    stop("All priors have to be named.")
  }
  unknown <- .unknown_priors(overrides, K)
  if (length(unknown) > 0) {
    stop(
      "Unknown prior(s): ", paste(unknown, collapse = ", "),
      ". Possible names are: ", paste(names(defaults), collapse = ", "), "."
    )
  }
  priors <- utils::modifyList(defaults, overrides)
  .check_priors(priors, K)
  priors
}

.default_priors <- function(K) {
  list(
    prior_Beta_loc = matrix(rep(0, K * K), nrow = K, ncol = K),
    prior_Beta_scale = matrix(rep(0.5, K * K), nrow = K, ncol = K),
    prior_Rho_loc = matrix(rep(0.5, K * K), nrow = K, ncol = K),
    prior_Rho_scale = matrix(rep(sqrt(0.5), K * K), nrow = K, ncol = K),
    prior_Rho_marginal = 0.25,
    prior_S = diag(K),
    prior_Eta = 1
  )
}

.unknown_priors <- function(priors, K) {
  setdiff(names(priors), names(.default_priors(K)))
}

# mirrors the constraints of the stan data blocks and prior densities
.check_priors <- function(priors, K) {
  mat_names <- c(
    "prior_Beta_loc", "prior_Beta_scale", "prior_Rho_loc", "prior_Rho_scale",
    "prior_S"
  )
  for (nm in mat_names) {
    x <- priors[[nm]]
    if (!is.matrix(x) || !is.numeric(x) || !all(dim(x) == K) || anyNA(x)) {
      stop("'", nm, "' must be a numeric ", K, " x ", K, " matrix.")
    }
  }
  for (nm in c("prior_Rho_marginal", "prior_Eta")) {
    x <- priors[[nm]]
    if (!is.numeric(x) || length(x) != 1 || is.na(x)) {
      stop("'", nm, "' must be a single number.")
    }
  }
  if (any(priors$prior_Beta_scale <= 0) || any(priors$prior_Rho_scale <= 0)) {
    stop("'prior_Beta_scale' and 'prior_Rho_scale' must be positive.")
  }
  if (any(priors$prior_Rho_loc <= 0 | priors$prior_Rho_loc >= 1)) {
    stop("'prior_Rho_loc' must be between 0 and 1.")
  }
  if (priors$prior_Rho_marginal <= 0 || priors$prior_Rho_marginal > 0.5) {
    stop("'prior_Rho_marginal' must be larger than 0 and at most 0.5.")
  }
  if (priors$prior_Eta < 1 || priors$prior_Eta != round(priors$prior_Eta)) {
    stop("'prior_Eta' must be a positive integer.")
  }
  s_ok <- isSymmetric(unname(priors$prior_S)) &&
    !is.null(tryCatch(chol(priors$prior_S), error = \(e) NULL))
  if (!s_ok) {
    stop("'prior_S' must be a symmetric positive definite matrix.")
  }
  invisible(priors)
}


#' Use the Posterior of a Fitted Model as Prior for the Temporal Effects
#'
#' @description Summarizes the posterior distribution of the temporal
#'   coefficients of a model fitted with \code{\link{stan_gvar}} by its means
#'   and standard deviations, in the format required for the \code{priors}
#'   argument of \code{\link{stan_gvar}}. This can be used to update a model
#'   with new data from the same person, or to inform the analysis of one
#'   data set with the results of another. Only priors for the temporal
#'   coefficients are returned; priors for the contemporaneous network are
#'   not supported.
#'
#' @param fit A \code{tsnet_fit} object obtained from
#'   \code{\link{stan_gvar}}.
#' @param inflate A positive number by which the posterior standard
#'   deviations are multiplied. Values larger than 1 make the resulting prior
#'   less informative. Default: \code{1}.
#' @param min_scale Lower bound for the resulting prior standard deviations.
#'   Default: \code{0.05}.
#'
#' @details The independent normal priors used by \code{\link{stan_gvar}}
#'   ignore posterior correlations between coefficients and deviations of the
#'   posterior from normality. The coefficients refer to standardized data
#'   (unless \code{center_only = TRUE} was used), which should be kept in mind
#'   if the variances of the variables differ between data sets.
#'
#' @return A list with the matrices \code{prior_Beta_loc} and
#'   \code{prior_Beta_scale}. Rows are outcomes and columns are lagged
#'   predictors, as expected by \code{\link{stan_gvar}}. Note that this is the
#'   transpose of the \code{beta_mu} matrix returned by
#'   \code{\link{stan_fit_convert}}.
#'
#' @seealso \code{\link{gvar_priors}}, \code{\link{stan_gvar}}
#'
#' @examples
#' \donttest{
#' data(ts_data)
#' fit_1 <- stan_gvar(ts_data[1:100, 1:3],
#'   cov_prior = "IW",
#'   n_chains = 2
#' )
#' new_priors <- priors_from_fit(fit_1, inflate = 2)
#' fit_2 <- stan_gvar(ts_data[101:200, 1:3],
#'   priors = new_priors,
#'   cov_prior = "IW",
#'   n_chains = 2
#' )
#' }
#'
#' @export
priors_from_fit <- function(fit, inflate = 1, min_scale = 0.05) {
  if (!inherits(fit, "tsnet_fit")) {
    stop("'fit' must be a 'tsnet_fit' object obtained from stan_gvar().")
  }
  if (!is.numeric(inflate) || length(inflate) != 1 || !isTRUE(inflate > 0)) {
    stop("'inflate' must be a single positive number.")
  }
  if (!is.numeric(min_scale) || length(min_scale) != 1 ||
    !isTRUE(min_scale > 0)) {
    stop("'min_scale' must be a single positive number.")
  }
  beta_draws <- .extract_draws(fit, "Beta")
  loc <- apply(beta_draws, c(2, 3), mean)
  scale <- pmax(apply(beta_draws, c(2, 3), stats::sd) * inflate, min_scale)
  cnames <- fit$arguments$cnames
  dimnames(loc) <- dimnames(scale) <- list(cnames, cnames)
  list(prior_Beta_loc = loc, prior_Beta_scale = scale)
}
