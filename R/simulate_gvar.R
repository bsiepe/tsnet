#' Simulate Data from a GVAR Model
#'
#' @description Simulates a multivariate time series from a lag-1 graphical
#'   vector autoregressive (GVAR) model with known parameters:
#'   \deqn{Y_t = \mu + B (Y_{t-1} - \mu) + \zeta_t, \quad
#'   \zeta_t \sim N(0, \Sigma).}
#'   The session random number generator is used, so call
#'   \code{\link{set.seed}} beforehand for reproducible results.
#'
#' @param Beta A \eqn{K \times K} matrix of temporal coefficients. Rows are
#'   outcomes and columns are lagged predictors, so that
#'   \code{Beta[i, j]} is the effect of variable \code{j} at time \eqn{t-1}
#'   on variable \code{i} at time \eqn{t}. This is the orientation of the
#'   \code{Beta} parameter in the Stan models and of \code{prior_Beta_loc} in
#'   \code{\link{stan_gvar}}. It is the transpose of the \code{beta_mu} matrix
#'   returned by \code{\link{stan_fit_convert}}.
#' @param Sigma A \eqn{K \times K} symmetric positive definite innovation
#'   covariance matrix.
#' @param n_t Number of time points to return.
#' @param burn_in Number of initial time points that are simulated and
#'   discarded, so that the returned series starts close to the stationary
#'   distribution. Default: \code{100}.
#' @param mu Optional vector of length \eqn{K} with the process means. The
#'   implied intercept is \eqn{(I - B)\mu}. Default: \code{NULL}, which
#'   corresponds to means of zero.
#'
#' @details The model has to be stationary, i.e., the spectral radius (largest
#'   absolute eigenvalue) of \code{Beta} has to be smaller than 1. Otherwise,
#'   an error is returned.
#'
#'   Note that \code{\link{stan_gvar}} standardizes the data before fitting.
#'   Parameter estimates are therefore comparable to the standardized
#'   parameters returned by \code{\link{gvar_stationary}}, not necessarily to
#'   \code{Beta} and \code{Sigma} themselves.
#'
#' @return A data frame with \code{n_t} rows and \eqn{K} columns that can be
#'   passed to \code{\link{stan_gvar}}. Column names are taken from
#'   \code{colnames(Beta)} if available and are \code{V1}, \code{V2}, ...
#'   otherwise. The data-generating parameters are attached as attribute
#'   \code{"truth"}, a list with elements \code{Beta}, \code{Sigma}, and
#'   \code{mu}.
#'
#' @seealso \code{\link{gvar_stationary}}, \code{\link{simulate.tsnet_fit}}
#'
#' @examples
#' Beta <- matrix(c(
#'   0.4, 0.2, 0.0,
#'   0.0, 0.3, 0.0,
#'   -0.2, 0.0, 0.5
#' ), nrow = 3, byrow = TRUE)
#' Sigma <- diag(3)
#' Sigma[1, 2] <- Sigma[2, 1] <- 0.3
#'
#' set.seed(1)
#' sim_data <- simulate_gvar(Beta, Sigma, n_t = 200)
#' head(sim_data)
#' attr(sim_data, "truth")$Beta
#'
#' @export
simulate_gvar <- function(Beta,
                          Sigma,
                          n_t,
                          burn_in = 100,
                          mu = NULL) {
  Sigma_chol <- .check_gvar_params(Beta, Sigma)
  K <- nrow(Beta)

  if (!.is_count(n_t) || n_t < 1) {
    stop("'n_t' must be a single positive integer.")
  }
  if (!.is_count(burn_in)) {
    stop("'burn_in' must be a single non-negative integer.")
  }
  if (is.null(mu)) {
    mu <- rep(0, K)
  }
  if (!is.numeric(mu) || length(mu) != K || !all(is.finite(mu))) {
    stop("'mu' must be a numeric vector with one element per variable.")
  }

  radius <- .spectral_radius(Beta)
  if (radius >= 1) {
    stop(
      "'Beta' does not imply a stationary process: its spectral radius is ",
      format(radius, digits = 3),
      ", but it has to be smaller than 1."
    )
  }

  Y <- .sim_var(Beta, Sigma_chol, n_t + burn_in)
  Y <- Y[seq_len(n_t) + burn_in, , drop = FALSE]
  Y <- sweep(Y, 2, mu, "+")

  cnames <- colnames(Beta)
  if (is.null(cnames)) {
    cnames <- paste0("V", seq_len(K))
  }
  out <- as.data.frame(Y)
  colnames(out) <- cnames
  attr(out, "truth") <- list(Beta = Beta, Sigma = Sigma, mu = mu)
  out
}

.is_count <- function(x) {
  is.numeric(x) && length(x) == 1 && !is.na(x) && x >= 0 && x == round(x)
}

.spectral_radius <- function(Beta) {
  max(Mod(eigen(Beta, only.values = TRUE)$values))
}

# returns the upper cholesky factor of Sigma so callers do not recompute it
.check_gvar_params <- function(Beta, Sigma) {
  if (!is.matrix(Beta) || !is.numeric(Beta) || nrow(Beta) != ncol(Beta)) {
    stop("'Beta' must be a square numeric matrix.")
  }
  K <- nrow(Beta)
  if (!is.matrix(Sigma) || !is.numeric(Sigma) || !all(dim(Sigma) == K)) {
    stop("'Sigma' must be a numeric matrix with the same dimensions as 'Beta'.")
  }
  if (!all(is.finite(Beta)) || !all(is.finite(Sigma))) {
    stop("'Beta' and 'Sigma' must not contain missing or infinite values.")
  }
  if (!isSymmetric(unname(Sigma))) {
    stop("'Sigma' must be symmetric.")
  }
  Sigma_chol <- tryCatch(chol(Sigma), error = \(e) NULL)
  if (is.null(Sigma_chol)) {
    stop("'Sigma' must be positive definite.")
  }
  Sigma_chol
}

# series starts at zero, without burn-in and without a stationarity check
.sim_var <- function(Beta, Sigma_chol, n) {
  K <- nrow(Beta)
  innov <- matrix(stats::rnorm(n * K), nrow = n, ncol = K) %*% Sigma_chol
  Y <- matrix(0, nrow = n, ncol = K)
  prev <- numeric(K)
  for (t in seq_len(n)) {
    prev <- drop(Beta %*% prev) + innov[t, ]
    Y[t, ] <- prev
  }
  Y
}


#' Stationary Distribution of a GVAR Model
#'
#' @description Computes the stationary (unconditional) covariance matrix
#'   implied by the parameters of a lag-1 GVAR model, together with the
#'   parameters of the same process after standardizing each variable to unit
#'   variance. As \code{\link{stan_gvar}} z-standardizes the data by default,
#'   the standardized parameters are the appropriate "truth" when checking how
#'   well a model recovers the parameters used in
#'   \code{\link{simulate_gvar}}.
#'
#' @inheritParams simulate_gvar
#'
#' @details The stationary covariance matrix \eqn{\Gamma_0} solves
#'   \eqn{\Gamma_0 = B \Gamma_0 B^\top + \Sigma}, which gives
#'   \eqn{vec(\Gamma_0) = (I - B \otimes B)^{-1} vec(\Sigma)}. With
#'   \eqn{D} the diagonal matrix of stationary standard deviations, the
#'   standardized process has temporal coefficients \eqn{D^{-1} B D} and
#'   innovation covariance \eqn{D^{-1} \Sigma D^{-1}}. Partial correlations
#'   of the innovations are not affected by the standardization.
#'
#' @return A list with elements
#'   \item{cov}{The stationary covariance matrix of the process.}
#'   \item{Beta_std}{Temporal coefficients of the standardized process
#'   (rows are outcomes, columns are lagged predictors).}
#'   \item{Sigma_std}{Innovation covariance matrix of the standardized
#'   process.}
#'   \item{pcor}{Partial correlations of the innovations (contemporaneous
#'   network), with a diagonal of zero.}
#'
#' @seealso \code{\link{simulate_gvar}}
#'
#' @examples
#' Beta <- matrix(c(
#'   0.4, 0.2, 0.0,
#'   0.0, 0.3, 0.0,
#'   -0.2, 0.0, 0.5
#' ), nrow = 3, byrow = TRUE)
#' Sigma <- diag(c(1, 2, 0.5))
#' Sigma[1, 2] <- Sigma[2, 1] <- 0.3
#'
#' gvar_stationary(Beta, Sigma)
#'
#' @export
gvar_stationary <- function(Beta, Sigma) {
  .check_gvar_params(Beta, Sigma)
  K <- nrow(Beta)
  radius <- .spectral_radius(Beta)
  if (radius >= 1) {
    stop(
      "'Beta' does not imply a stationary process: its spectral radius is ",
      format(radius, digits = 3),
      ", but it has to be smaller than 1."
    )
  }

  G0 <- solve(diag(K^2) - Beta %x% Beta, as.vector(Sigma)) |>
    matrix(nrow = K, ncol = K)
  G0 <- (G0 + t(G0)) / 2
  sds <- sqrt(diag(G0))

  Theta <- solve(Sigma)
  pcor <- -stats::cov2cor(Theta)
  diag(pcor) <- 0

  out <- list(
    cov = G0,
    Beta_std = Beta * outer(1 / sds, sds),
    Sigma_std = Sigma / outer(sds, sds),
    pcor = pcor
  )
  lapply(out, \(m) `dimnames<-`(m, dimnames(Beta)))
}
