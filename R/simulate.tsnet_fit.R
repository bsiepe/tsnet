#' Simulate Data from a Fitted GVAR Model
#'
#' @description Simulates new time series from the draws of a model fitted
#'   with \code{\link{stan_gvar}}. Each simulated data set is based on a
#'   different draw of the temporal coefficients and the innovation covariance
#'   matrix and is generated with \code{\link{simulate_gvar}}. For a regular
#'   fit, this results in posterior predictive simulations. For a fit obtained
#'   with \code{prior_only = TRUE}, this results in prior predictive
#'   simulations.
#'
#' @param object A \code{tsnet_fit} object obtained from
#'   \code{\link{stan_gvar}}.
#' @param nsim Number of data sets to simulate. Each data set uses a different
#'   draw, which is selected at random without replacement. Default: \code{1}.
#' @param seed Either \code{NULL} (default) or an integer that is used in a
#'   call to \code{\link{set.seed}} before simulating. If set, the state of
#'   the random number generator is restored afterwards, as in other
#'   \code{\link[stats]{simulate}} methods.
#' @param n_t Number of time points per simulated data set. Defaults to the
#'   number of time points of the data used for fitting.
#' @param nonstationary How to treat draws with non-stationary temporal
#'   coefficients (spectral radius of at least 1). \code{"drop"} (default)
#'   only uses stationary draws, so that the simulated data are conditional on
#'   stationarity; a message reports how many draws were not used.
#'   \code{"keep_flagged"} uses all draws and
#'   flags the non-stationary ones in the \code{"draws"} attribute; data for
#'   these draws are simulated from a starting value of zero without burn-in
#'   and will typically diverge. \code{"error"} stops if one of the selected
#'   draws is non-stationary.
#' @param ... Further arguments passed to \code{\link{simulate_gvar}} for
#'   stationary draws, currently only \code{burn_in}.
#'
#' @details The simulated data are on the scale that was used for fitting:
#'   \code{\link{stan_gvar}} z-standardizes the data (or only centers them if
#'   \code{center_only = TRUE}) and estimates a model without intercepts, and
#'   the original means and standard deviations are not stored. Simulated
#'   series therefore have process means of zero and are not transformed back
#'   to the original scale.
#'
#'   If the model was fitted with \code{rmv_overnight = TRUE}, the day
#'   structure is ignored and a single continuous time series is simulated.
#'
#'   Non-stationary draws can be frequent, in particular for prior-only fits
#'   with the default priors and for models with many variables. Draws are
#'   never reused. If fewer than \code{nsim} usable draws are available, fewer
#'   data sets are returned with a warning. Draws with non-finite values or
#'   a numerically non-positive-definite covariance matrix are never used.
#'
#' @return A list of data frames, one per simulated data set, as returned by
#'   \code{\link{simulate_gvar}}. The parameter draw used for each data set is
#'   attached to it as attribute \code{"truth"}. The list has the following
#'   attributes:
#'   \item{seed}{The state of the random number generator, as in other
#'   \code{\link[stats]{simulate}} methods.}
#'   \item{draws}{A data frame with the index of the draw used for each data
#'   set, its spectral radius, and whether it is stationary. The index refers
#'   to the (permuted) order of draws returned by
#'   \code{\link[rstan:extract]{rstan::extract}}.}
#'   \item{n_draws}{A named vector with the number of draws in the fit
#'   (\code{total}), the number of \code{stationary} and
#'   \code{nonstationary} draws, the number of \code{invalid} draws
#'   (non-finite values or a covariance matrix that is not positive
#'   definite), and the number of \code{requested} and \code{returned}
#'   data sets.}
#'
#' @seealso \code{\link{simulate_gvar}}, \code{\link{stan_gvar}}
#'
#' @examples
#' \donttest{
#' data(ts_data)
#' example_data <- ts_data[1:100, 1:3]
#' fit <- stan_gvar(example_data,
#'   cov_prior = "IW",
#'   n_chains = 2
#' )
#' sims <- simulate(fit, nsim = 5, seed = 1)
#' head(sims[[1]])
#' attr(sims, "n_draws")
#' }
#'
#' @importFrom stats simulate
#' @export
simulate.tsnet_fit <- function(object,
                               nsim = 1,
                               seed = NULL,
                               n_t = NULL,
                               nonstationary = c(
                                 "drop", "keep_flagged", "error"
                               ),
                               ...) {
  nonstationary <- match.arg(nonstationary)
  if (!.is_count(nsim) || nsim < 1) {
    stop("'nsim' must be a single positive integer.")
  }
  if (is.null(n_t)) {
    n_t <- object$arguments$n_t
  }
  if (!.is_count(n_t) || n_t < 1) {
    stop("'n_t' must be a single positive integer.")
  }

  # same seed handling as stats:::simulate.lm
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    stats::runif(1)
  }
  if (is.null(seed)) {
    rng_state <- get(".Random.seed", envir = .GlobalEnv)
  } else {
    old_seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    rng_state <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv))
  }

  beta_draws <- .extract_draws(object, "Beta")
  sigma_draws <- .extract_draws(object, "Sigma")
  n_total <- dim(beta_draws)[1]

  radius <- .draw_radius(beta_draws)
  sigma_chol <- lapply(seq_len(n_total), function(s) {
    Sigma <- .draw_matrix(sigma_draws, s)
    if (!all(is.finite(Sigma))) {
      return(NULL)
    }
    tryCatch(chol(.symmetrize(Sigma)), error = \(e) NULL)
  })
  valid <- !is.na(radius) & !vapply(sigma_chol, is.null, logical(1))
  stationary <- valid & radius < 1

  n_nonstationary <- sum(valid & !stationary)
  if (nonstationary == "drop" && n_nonstationary > 0) {
    message(
      n_nonstationary, " of ", sum(valid), " draws are non-stationary and ",
      "are not used, so the simulated data are conditional on stationarity. ",
      "Use nonstationary = \"keep_flagged\" to include them."
    )
  }

  usable <- which(valid & (stationary | nonstationary != "drop"))
  n_use <- min(nsim, length(usable))
  ids <- usable[sample.int(length(usable), n_use)]

  if (nonstationary == "error" && !all(stationary[ids])) {
    stop(
      sum(!stationary[ids]), " of the ", n_use,
      " selected draws are non-stationary."
    )
  }
  if (n_use < nsim) {
    warning(
      "Only ", n_use, " usable draws are available, so fewer than the ",
      "requested ", nsim, " data sets are returned."
    )
  }

  cnames <- object$arguments$cnames
  sims <- lapply(ids, function(s) {
    Beta <- .draw_matrix(beta_draws, s)
    Sigma <- .symmetrize(.draw_matrix(sigma_draws, s))
    dimnames(Beta) <- dimnames(Sigma) <- list(cnames, cnames)
    if (stationary[s]) {
      return(simulate_gvar(Beta, Sigma, n_t = n_t, ...))
    }
    Y <- as.data.frame(.sim_var(Beta, sigma_chol[[s]], n_t))
    if (!is.null(cnames)) {
      colnames(Y) <- cnames
    }
    attr(Y, "truth") <- list(Beta = Beta, Sigma = Sigma, mu = rep(0, ncol(Y)))
    Y
  })

  attr(sims, "seed") <- rng_state
  attr(sims, "draws") <- data.frame(
    draw = ids,
    spectral_radius = radius[ids],
    stationary = stationary[ids]
  )
  attr(sims, "n_draws") <- c(
    total = n_total,
    stationary = sum(stationary),
    nonstationary = n_nonstationary,
    invalid = sum(!valid),
    requested = nsim,
    returned = n_use
  )
  sims
}

# draws in stan orientation: [draw, row, column] with rows of Beta as outcomes
.extract_draws <- function(fit, par) {
  rstan::extract(fit$stan_fit, pars = par)[[par]]
}

# NA for draws with non-finite values
.draw_radius <- function(beta_draws) {
  vapply(seq_len(dim(beta_draws)[1]), function(s) {
    Beta <- .draw_matrix(beta_draws, s)
    if (all(is.finite(Beta))) .spectral_radius(Beta) else NA_real_
  }, numeric(1))
}

# keeps the matrix structure for models with a single variable
.draw_matrix <- function(draws, s) {
  K <- dim(draws)[2]
  matrix(draws[s, , ], nrow = K, ncol = K)
}

# stan draws of a covariance matrix are only symmetric up to rounding
.symmetrize <- function(m) {
  (m + t(m)) / 2
}
