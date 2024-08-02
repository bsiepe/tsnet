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
#'   covariance matrix. Options are "LKJ" or "IW" (Inverse-Wishart). Default: \code{"LKJ"}.
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
#'   scale) the data. Default is \code{FALSE}.
#' @param ... Additional arguments passed to the \code{\link[rstan::sampling]{rstan::sampling}} or
#'   \code{\link[rstan::vb]{rstan::vb}}  function.
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
#'   interval of \eqn{[0,1]} and then assigning a proportional (mean-variance
#'   parameterized) beta prior:
#'   \deqn{(\rho_{ij}+1)/2 \sim Beta_{prop}(PriorRhoLoc, PriorRhoScale)}
#'   A beta location parameter of 0.5 translates to an expected correlation of 0.
#'   The variance parameter of \eqn{\sqrt(0.5)} implies a uniform distribution of
#'   partial correlations.
#'   The user can specify a different prior distribution by a matrix
#'   \code{prior_Rho_loc} and a matrix \code{prior_Rho_scale} with the same dimensions as
#'   the partial correlation matrix. Additionally, the user can change \eqn{eta}
#'   via the \code{prior_Eta} parameter.
#'
#'   The Inverse-Wishart prior is a distribution on the innovation covariance
#'   matrix \eqn{\Sigma}:
#'   \deqn{\Sigma \sim IW(\nu, S)}
#'   where \eqn{\nu} is the degrees of freedom and \eqn{S} is the scale matrix. We here
#'   use the default prior of \deqn{nu = delta + p - 1} for the degrees of freedom,
#'   where \eqn{\delta} is defined as \eqn{s_{\rho}^{-1}-1} and \eqn{s_{\rho}} is the
#'   standard deviation of the implied marginal beta distribution of the
#'   partial correlations. For the scale matrix \eqn{S}, we use the identity matrix
#'   \eqn{I_p} of order \eqn{p}.
#'   The user can set a prior on the expected standard deviation of the partial
#'   correlations by specifying a \code{prior_Rho_marginal} parameter. The default
#'   value is 0.25, which has worked well in a simulation study.
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
#'  \item{arguments}{The number of variables "p", the number of time points "n_t", the column names "cnames", and the arguments used in the function call.}
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

stan_gvar <-
  function(data,
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
           ...) {
    if(isTRUE(center_only)){
      Y <- apply(data, MARGIN = 2, scale, center = TRUE, scale = FALSE)
    } else{
      Y <- apply(data, MARGIN = 2, scale, center = TRUE, scale = TRUE)
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

    if(is.null(beep)){
      beep <- seq(1, n_t)
    }


    # Specify Priors
    if(is.null(priors[["prior_Beta_loc"]])){
        prior_Beta_loc <- matrix(rep(0, K*K), nrow = K, ncol = K)
    } else {
        prior_Beta_loc <- priors[["prior_Beta_loc"]]
    }

    if(is.null(priors[["prior_Beta_scale"]])){
        prior_Beta_scale <- matrix(rep(.5, K*K), nrow = K, ncol = K)
    } else {
        prior_Beta_scale <- priors[["prior_Beta_scale"]]
    }

    if(is.null(priors[["prior_Rho_loc"]])){
        prior_Rho_loc <- matrix(rep(.5, K*K), nrow = K, ncol = K)
    } else {
        prior_Rho_loc <- priors[["prior_Rho_loc"]]
    }

    if(is.null(priors[["prior_Rho_scale"]])){
        prior_Rho_scale <- matrix(rep(sqrt(.5), K*K), nrow = K, ncol = K)
    } else {
        prior_Rho_scale <- priors[["prior_Rho_scale"]]
    }

    if(is.null(priors[["prior_Rho_marginal"]])){
        prior_Rho_marginal <- 0.25
    } else {
        prior_Rho_marginal <- priors[["prior_Rho_marginal"]]
    }

    if(is.null(priors[["prior_S"]])){
      prior_S <- diag(K)
    } else {
      prior_S <- priors[["prior_S"]]
    }

    if(is.null(priors[["prior_Eta"]])){
      prior_Eta <- 1
    } else {
      prior_Eta <- priors[["prior_Eta"]]
    }

    # Convert SD to delta: SD = 1/(delta+1), delta = (1 / SD) - 1
    prior_delta <- (1 / prior_Rho_marginal) - 1



    # Choose model to fit
    if (cov_prior == "LKJ") {
      # Stan Data
      stan_data <- list(
        K = K,
        "T" = n_t,
        Y = as.matrix(Y),
        beep = beep,
        prior_Rho_loc = prior_Rho_loc,
        prior_Rho_scale = prior_Rho_scale,
        prior_Beta_loc = prior_Beta_loc,
        prior_Beta_scale = prior_Beta_scale,
        prior_Eta = prior_Eta
      )

      if (isTRUE(rmv_overnight)) {
        # remove overnight effects
        model_name <- "VAR_LKJ_beep"
      } else{
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
        prior_Beta_loc = prior_Beta_loc,
        prior_Beta_scale = prior_Beta_scale,
        prior_S = prior_S,
        prior_delta = prior_delta
      )

      if (isTRUE(rmv_overnight)) {
        # remove overnight effects
        model_name <- "VAR_wishart_beep"
      } else{
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

        # Run sampler
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
                 fn_args = all_args
                 )

    ret_fit <- list(
      stan_fit = stan_fit,
      arguments = args
    )

    class(ret_fit) <- c("tsnet_fit", class(ret_fit))
    return(ret_fit)

}
