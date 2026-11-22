#' Fit Bayesian Graphical Vector Autoregressive (GVAR) Models with Stan
#'
#' @description This function fits a Bayesian GVAR model to the provided data
#'   using Stan. The estimation procedure is described further in Siepe, Kloft &
#'   Heck (2023) <doi:10.31234/osf.io/uwfjc>. The current implementation allows
#'   for a normal prior on the temporal effects and either an Inverse Wishart or
#'   an LKJ prior on the contemporaneous effects. \code{rstan} is used as a backend
#'   for fitting the model in Stan. Data should be provided in long format, where
#'   the columns represent the variables and the rows represent the time points.
#'   Data are automatically z-scaled for estimation. The model currently does not
#'   support missing data.
#'
#' @param data A data frame or matrix containing the time series data of a
#'   single subject. The data should be in long format,  where the columns
#'   represent the variables and the rows represent the time points. See the
#'   example data \code{\link{ts_data}} for the correct format.
#' @param beep A vector of beeps with length of \code{nrow(data)}. The beep indicator
#'   can be used to remove overnight effects from the last beep of a day to the
#'   first beep of the next day. This should be a vector of positive integers.
#'   If left empty, the function will assume that there are no overnight
#'   effects to remove.
#' @param priors A list of prior distributions for the model parameters. This
#'   should be a named list, with names corresponding to the parameter names and
#'   values corresponding to the prior distributions. The following priors can
#'   be specified:
#'   \itemize{
#'   \item \code{prior_Beta_loc} A matrix of the same dimensions as the beta matrix
#'   `B` containing the mean of the prior distribution for the beta coefficients.
#'   \item \code{prior_Beta_scale} A matrix of the same dimensions as the beta matrix
#'   `B` containing the standard deviation of the prior distribution for the beta
#'   coefficients.}
#'
#' @param method A string indicating the method to use for fitting the model.
#'   Options are "sampling" (for MCMC estimation) or "variational" (for
#'   variational inference). We currently recommend only using MCMC estimation.
#'   Default: \code{"sampling"}.
#' @param cov_prior A string indicating the prior distribution to use for the
#'   covariance matrix. Options are "LKJ" or "IW" (Inverse-Wishart). Default: \code{"IW"}.
#' @param rmv_overnight A logical indicating whether to remove overnight
#'   effects. Default is \code{FALSE}. If `TRUE`, the function will remove overnight
#'   effects from the last beep of a day to the first beep of the next day.
#'   This requires the \code{beep} argument to be specified.
#' @param iter_sampling An integer specifying the number of iterations for the
#'   sampling method. Default is \code{500}.
#' @param iter_warmup An integer specifying the number of warmup iterations for
#'   the sampling method. Default is \code{500}.
#' @param n_chains An integer specifying the number of chains for the sampling
#'   method. Default is 4. If variational inference is used, the number of
#'   iterations is calculated as \code{iter_sampling}*\code{n_chains}.
#' @param n_cores An integer specifying the number of cores to use for parallel
#'   computation. Default is \code{1}. \code{rstan} is used for parallel computation.
#' @param center_only A logical indicating whether to only center (and not
#'   scale) the data. Default is \code{FALSE}, which z-standardizes each
#'   variable with its sample mean and standard deviation, so that the
#'   temporal coefficients are standardized coefficients. The default priors
#'   are specified for standardized data.
#' @param return_all A logical indicating whether to return all model inputs, including
#' the data and prior objects. Default is \code{TRUE}.
#' @param ahead An integer specifying the forecast horizon. Default is \code{0}. If `ahead` is
#'  greater than 0, the function will return posterior predictive forecasts for the specified number of time points.
#'  This functionality is experimental and may not work as expected.
#' @param ... Additional arguments passed to the \code{\link[rstan:sampling]{rstan::sampling}} or
#'   \code{\link[rstan:vb]{rstan::vb}}  function.
#' @param prior_only A logical indicating whether to sample from the prior
#'   distribution only, ignoring the likelihood of the data. Default is
#'   \code{FALSE}. The data are still required to determine the number and
#'   names of the variables and the number of time points, but their values
#'   are not used, so that placeholder data such as
#'   \code{matrix(0, 100, 3)} can be passed. The resulting draws can be used
#'   to inspect the
#'   priors (e.g., with \code{\link{posterior_plot}}) and, together with
#'   \code{\link{simulate.tsnet_fit}}, for prior predictive simulation. Only
#'   available for \code{method = "sampling"} and \code{ahead = 0}.
#'
#' @details \bold{General Information}
#'
#'   In a Graphical Vector Autoregressive (GVAR) model of lag 1, each variable
#'   is regressed on itself and all other variables at the previous timepoint to
#'   obtain estimates of the temporal association between variables
#'   (encapsulated in the beta matrix). This is the "Vector Autoregressive" part
#'   of the model. Additionally, the innovation structure at each time point
#'   (which resembles the residuals) is modeled to obtain estimates of the
#'   contemporaneous associations between all variables (controlling for the
#'   lagged effects). This is typically represented in the partial correlation
#'   (pcor) matrix. If the model is represented and interpreted as a network,
#'   variables are called
#' \emph{nodes}, \emph{edges} represent the statistical association between the nodes, and
#' \emph{edge weights} quantify the strength of these associations.
#'
#'   \bold{Model}
#'
#'   Let \eqn{Y} be a matrix with \eqn{n} rows and \eqn{p} columns, where \eqn{n_t} is the
#'   number of time points and \eqn{p} is the number of variables. The GVAR model is
#'   given by the following equations: \deqn{Y_t = B* Y_{t-1} + \zeta_t}
#'   \deqn{\zeta_t \sim  N(0, \Sigma)} where \eqn{B} is a \eqn{p \times p} matrix of VAR
#'   coefficients between variables \eqn{i} and \eqn{j} (\eqn{\beta_{ij}}), \eqn{\zeta_t} contains the
#'   innovations at time point \eqn{t}, and \eqn{\Sigma} is a \eqn{p \times p} covariance matrix.
#'   The inverse of \eqn{\Sigma} is the precision matrix, which is used to obtain the
#'   partial correlations between variables (\eqn{\rho_{ij}}). The model setup is
#'   explained in more detail in Siepe, Kloft & Heck (2023)
#'   <doi:10.31234/osf.io/uwfjc>.
#'
#'
#'   \bold{Prior Setup}
#'
#'   For the \eqn{p \times p} temporal matrix B (containing the \eqn{\beta} coefficients), we use
#'   a normal prior distribution on each individual parameter: \deqn{\beta_{ij}
#'   \sim N(PriorBetaLoc_{ij}, PriorBetaScale_{ij})} where \code{PriorBetaLoc} is the
#'   mean of the prior distribution and \code{PriorBetaScale} is the standard
#'   deviation of the prior distribution. The default prior is a weakly
#'   informative normal distribution with mean 0 and standard deviation 0.5. The
#'   user can specify a different prior distribution by a matrix
#'   \code{prior_Beta_loc} and a matrix \code{prior_Beta_scale} with the same dimensions
#'   as \eqn{B}.
#'
#'
#'   Both a Lewandowski-Kurowicka-Joe (LKJ) and an Inverse-Wishart (IW)
#'   distribution can be used as a prior for the contemporaneous network.
#'   However, the LKJ prior does not allow for direct specifications of priors
#'   on the partial correlations. We implemented a workaround to enable priors
#'   on specific partial correlations (described below). We consider this
#'   feature experimental would advise users wishing to implement edge-specific
#'   priors in the contemporaneous network to preferentially use IW priors.
#'
#'   The LKJ prior is a distribution on the correlation matrix, which is
#'   parameterized by the shape parameter \eqn{\eta}. To enable edge-specific priors
#'   on the partial correlations, we use the workaround of a "joint" prior
#'   that, in addition to the LKJ on the correlation matrix itself, allows for
#'   an additional beta prior on each of the partial correlations. We first
#'   assigned an uninformed LKJ prior to the Cholesky factor decomposition of
#'   the correlation matrix of innovations: \deqn{\Omega_L \sim
#'   LKJ-Cholesky(\eta)}
#'   For \eqn{\eta = 1}, this implies a symmetric marginal
#'   scaled beta distribution on the zero-order correlations \eqn{\omega_{ij}}.
#'   \deqn{(\omega_{ij}+1)/2 \sim Beta(p/2, p/2)}
#'   We can then obtain the covariance matrix and,
#'   subsequently, the precision matrix (see Siepe, Kloft & Heck, 2023,
#'   for details).
#'   The second part of the prior is a beta prior on each partial correlation
#'   \eqn{\rho_{ij}} (obtained from the off-diagonal elements of the precision matrix).
#'   This prior was assigned by transforming the partial correlations to the
#'   interval of \eqn{[0,1]} and then assigning a beta prior in the
#'   mean-precision parameterization (\code{beta_proportion} in Stan):
#'   \deqn{(\rho_{ij}+1)/2 \sim Beta(\mu_{ij} \kappa_{ij},
#'   (1 - \mu_{ij}) \kappa_{ij})}
#'   where \eqn{\mu_{ij}} is given by \code{prior_Rho_loc} and
#'   \eqn{\kappa_{ij}} by \code{prior_Rho_scale}. A location of 0.5 centers
#'   the beta prior at a partial correlation of 0. Note that
#'   \code{prior_Rho_scale} is a precision and not a variance or standard
#'   deviation: larger values lead to a more concentrated prior. A value of
#'   2 (with a location of 0.5), which is the default, corresponds to a
#'   uniform beta prior. In this case the partial correlations marginally
#'   follow \eqn{(\rho_{ij}+1)/2 \sim Beta(\eta, \eta)}, which is uniform
#'   for the default of \eqn{\eta = 1}. Beta priors with a shape parameter
#'   below 1 (i.e., \eqn{\mu_{ij} \kappa_{ij} < 1} or
#'   \eqn{(1 - \mu_{ij}) \kappa_{ij} < 1}) are U-shaped and unbounded at
#'   partial correlations of -1 or 1. They favor very strong partial
#'   correlations, and sampling from such a prior without data
#'   (\code{prior_only = TRUE}) is unreliable. Note that versions up to
#'   0.2.0 used a default of \eqn{\sqrt{0.5}}, which is such a prior
#'   (\eqn{Beta(0.354, 0.354)}); set \code{prior_Rho_scale} to this value
#'   to reproduce earlier results.
#'   The user can specify a different prior distribution by a matrix
#'   \code{prior_Rho_loc} and a matrix \code{prior_Rho_scale} with the same dimensions as
#'   the partial correlation matrix, of which only the upper triangle is
#'   used. Additionally, the user can change \eqn{\eta}
#'   via the \code{prior_Eta} parameter, which has to be a positive integer.
#'
#'   With \code{cov_prior = "IW"}, an Inverse-Wishart prior is placed on the
#'   precision matrix of the innovations, \eqn{\Theta = \Sigma^{-1}}:
#'   \deqn{\Theta \sim IW(\nu, S)}
#'   where \eqn{\nu} is the degrees of freedom and \eqn{S} is the scale matrix.
#'   This is equivalent to a Wishart prior on the innovation covariance
#'   matrix, \eqn{\Sigma \sim W(\nu, S^{-1})}. The degrees of freedom are
#'   \deqn{\nu = \delta + p - 1}
#'   where \eqn{\delta} is defined as \eqn{v_{\rho}^{-1}-1} and \eqn{v_{\rho}}
#'   is the value of \code{prior_Rho_marginal}. For the scale matrix \eqn{S},
#'   we use the identity matrix \eqn{I_p} of order \eqn{p} by default. In this
#'   case, the partial correlations marginally follow
#'   \eqn{(\rho_{ij}+1)/2 \sim Beta(\delta/2, \delta/2)}, which has a variance
#'   of \eqn{1/(\delta+1) = v_{\rho}} on the scale of the partial
#'   correlations. \code{prior_Rho_marginal} is therefore the prior variance
#'   of the partial correlations, and its square root is their prior standard
#'   deviation. The default value of 0.25 (\eqn{\delta = 3}) corresponds to a
#'   prior standard deviation of 0.5. Values have to be larger than 0 and at
#'   most 0.5.
#'   Additionally, the user can specify a \code{prior_S} parameter to set a different
#'   scale matrix.
#'
#'   \bold{Sampling}
#'   The model can be fitted using either MCMC sampling or variational
#'   inference via \code{rstan}. Per default, the model is fitted using the Stan
#'   Hamiltonian Monte Carlo (HMC) No U-Turn (NUTS) sampler with 4 chains,
#'   500 warmup iterations and 500 sampling iterations. We use a default
#'   target average acceptance probability \code{adapt_delta} of 0.8. As the output
#'   is returned as a standard \code{stanfit} object, the user can use the
#'   \code{rstan} package to extract and analyze the results and obtain convergence
#'   diagnostics.
#'
#'
#'
#' @return A \code{tsnet_fit} object in list format. The object contains the
#'  following elements:
#'  \item{fit}{A stanfit object containing the fitted model.}
#'  \item{arguments}{The number of variables "p", the number of time points "n_t", the column names "cnames", and the arguments used in the function call. "center" and "scale" contain the sample means and standard deviations used to standardize the data (standard deviations of 1 if \code{center_only = TRUE}, both \code{NULL} for prior-only fits).}
#'
#' @importFrom rstan sampling vb
#' @importFrom utils modifyList
#' @examples
#' \donttest{
#' # Load example data
#' data(ts_data)
#' example_data <- ts_data[1:100,1:3]
#'
#' # Fit the model
#' fit <- stan_gvar(example_data,
#'                  method = "sampling",
#'                  cov_prior = "IW",
#'                  n_chains = 2)
#' print(fit)
#'}
#' @export

stan_gvar <- function(data,
                      beep = NULL,
                      priors = NULL,
                      method = "sampling",
                      cov_prior = "IW",
                      rmv_overnight = FALSE,
                      iter_sampling = 500,
                      iter_warmup = 500,
                      n_chains = 4,
                      n_cores = 1,
                      center_only = FALSE,
                      return_all = TRUE,
                      ahead = 0,
                      ...,
                      prior_only = FALSE) {
  if (!is.logical(prior_only) || length(prior_only) != 1 ||
    is.na(prior_only)) {
    stop("'prior_only' must be either TRUE or FALSE.")
  }
  if (prior_only && ahead > 0) {
    stop("'prior_only = TRUE' cannot be combined with 'ahead > 0'.")
  }
  if (prior_only && method == "variational") {
    stop("'prior_only = TRUE' is only available for method = \"sampling\".")
  }

  if (prior_only) {
    # the likelihood is skipped, so only the dimensions of the data matter.
    # this also allows constant placeholder data
    Y <- matrix(0, nrow(data), ncol(data), dimnames = dimnames(data))
    data_center <- NULL
    data_scale <- NULL
  } else if (isTRUE(center_only)) {
    Y <- apply(data, MARGIN = 2, scale, center = TRUE, scale = FALSE)
    data_center <- colMeans(data)
    data_scale <- rep(1, ncol(data))
  } else {
    Y <- apply(data, MARGIN = 2, scale, center = TRUE, scale = TRUE)
    data_center <- colMeans(data)
    data_scale <- apply(data, MARGIN = 2, stats::sd)
  }

  K <- ncol(data)
  n_t <- nrow(data)
  cnames <- colnames(data)

  # Store arguments
  mc <- match.call()
  fl <- formals()
  missing_args <- setdiff(names(fl), names(mc))
  missing_args <- Map(function(arg, default)
    if (!is.null(default)) mc[[arg]] <- default, missing_args, fl[missing_args])
  all_args <- c(as.list(mc), missing_args)

  if (is.null(beep)) {
    beep <- seq(1, n_t)
  }

  # Initialize priors with default values
  default_priors <- .default_priors(K)

  # Update default priors with any custom priors provided by the user
  if (!is.null(priors)) {
    if (is.null(names(priors)) || any(names(priors) == "")) {
      warning("'priors' must be a named list. Unnamed entries are ignored.")
    }
    unknown_priors <- .unknown_priors(priors, K)
    if (length(unknown_priors) > 0) {
      warning(
        "Unknown prior(s) will be ignored: ",
        paste(unknown_priors, collapse = ", "),
        ". Possible names are: ",
        paste(names(default_priors), collapse = ", "), "."
      )
    }
    priors <- modifyList(default_priors, priors)
  } else {
    priors <- default_priors
  }

  # prior_Rho_marginal acts as the prior variance of the partial correlations
  # (for prior_S = I): var = 1 / (delta + 1), delta = (1 / var) - 1
  prior_delta <- (1 / priors$prior_Rho_marginal) - 1

  # Forecasting
  if (ahead > 0){
    compute_log_lik <- 1
  } else {
    compute_log_lik <- 0
  }

  # Choose model to fit
  if (cov_prior == "LKJ") {
    if (prior_only && .unbounded_rho_prior(priors)) {
      warning(
        "Prior-only sampling with 'cov_prior = \"LKJ\"' is unreliable if the ",
        "beta prior on a partial correlation is unbounded at -1 or 1. ",
        "Consider values of at least 2 for ",
        "'prior_Rho_scale' (with 'prior_Rho_loc' of 0.5)."
      )
    }
    # Stan Data
    stan_data <- list(
      K = K,
      "T" = n_t,
      Y = as.matrix(Y),
      beep = beep,
      prior_Rho_loc = priors$prior_Rho_loc,
      prior_Rho_scale = priors$prior_Rho_scale,
      prior_Beta_loc = priors$prior_Beta_loc,
      prior_Beta_scale = priors$prior_Beta_scale,
      prior_Eta = priors$prior_Eta,
      ahead = ahead,
      Y_future = matrix(rep(0, K), ncol = K),
      compute_log_lik = compute_log_lik,
      prior_only = as.integer(prior_only)
    )

    if (isTRUE(rmv_overnight)) {
      # remove overnight effects
      model_name <- "VAR_LKJ_beep"
    } else {
      # standard model
      model_name <- "VAR_LKJ"
    }
  }
  if (cov_prior == "IW") {
    # Stan Data
    stan_data <- list(
      K = K,
      "T" = n_t,
      Y = as.matrix(Y),
      beep = beep,
      prior_Beta_loc = priors$prior_Beta_loc,
      prior_Beta_scale = priors$prior_Beta_scale,
      prior_S = priors$prior_S,
      prior_delta = prior_delta,
      ahead = ahead,
      Y_future = matrix(rep(0, K), ncol = K),
      compute_log_lik = compute_log_lik,
      prior_only = as.integer(prior_only)
    )

    if (isTRUE(rmv_overnight)) {
      # remove overnight effects
      model_name <- "VAR_wishart_beep"
    } else {
      # standard model
      model_name <- "VAR_wishart"
    }
  }

  if (method == "sampling") {
    default_args <- list(
      object = stanmodels[[model_name]],
      data = stan_data,
      chains = n_chains,
      cores = n_cores,
      iter = iter_sampling + iter_warmup,
      warmup = iter_warmup,
      refresh = 500,
      thin = 1,
      init = .1,
      control = list(adapt_delta = .8)
    )

    stan_fit <- do.call(rstan::sampling,
                        utils::modifyList(default_args, list(...)))

  }
  if (method == "variational") {
    default_args <- list(
      object = stanmodels[[model_name]],
      data = stan_data,
      init = .1,
      tol_rel_obj = .001,
      output_samples = iter_sampling * n_chains
    )

    stan_fit <- do.call(rstan::vb,
                        utils::modifyList(default_args, list(...)))
  }


    args <- list(p = K,
                 n_t = n_t,
                 cnames = cnames,
                 fn_args = all_args,
                 prior_only = prior_only,
                 center = data_center,
                 scale = data_scale
                 )
    if(isTRUE(return_all)){
      args$data <- stan_data$Y
      args$priors <- priors
    }


    ret_fit <- list(
      stan_fit = stan_fit,
      arguments = args
    )

    class(ret_fit) <- c("tsnet_fit", class(ret_fit))
    return(ret_fit)

}
