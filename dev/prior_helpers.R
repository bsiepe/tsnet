# helpers shared by dev/prior_checks.R and dev/prior_spec_check.R

probs <- c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95)

# LKJ(eta) draws via the c-vine method of Lewandowski, Kurowicka & Joe (2009)
# returns an array [draw, K, K] of correlation matrices
rlkj <- function(n, K, eta) {
  P <- array(0, c(n, K, K))
  S <- array(0, c(n, K, K))
  for (i in seq_len(K)) S[, i, i] <- 1
  b <- eta + (K - 1) / 2
  for (k in seq_len(K - 1)) {
    b <- b - 0.5
    for (i in (k + 1):K) {
      P[, k, i] <- 2 * rbeta(n, b, b) - 1
      p <- P[, k, i]
      if (k > 1) {
        for (l in (k - 1):1) {
          p <- p * sqrt((1 - P[, l, i]^2) * (1 - P[, l, k]^2)) +
            P[, l, i] * P[, l, k]
        }
      }
      S[, k, i] <- S[, i, k] <- p
    }
  }
  S
}

# partial correlations for an array [draw, K, K] of covariance or
# precision matrices
pcor_draws <- function(x, is_precision = FALSE) {
  out <- array(0, dim(x))
  for (s in seq_len(dim(x)[1])) {
    Theta <- if (is_precision) {
      x[s, , ]
    } else {
      tryCatch(solve(x[s, , ]), error = function(e) NULL)
    }
    # singular matrices, e.g. from chains stuck at a correlation of 1
    if (is.null(Theta)) {
      out[s, , ] <- NA
      next
    }
    d <- 1 / sqrt(diag(Theta))
    R <- -Theta * outer(d, d)
    diag(R) <- 0
    out[s, , ] <- R
  }
  out
}

kish_ess <- function(w) sum(w)^2 / sum(w^2)

wquantile <- function(x, w, probs) {
  o <- order(x)
  cw <- cumsum(w[o]) / sum(w)
  x[o][vapply(probs, \(p) which(cw >= p)[1], integer(1))]
}

# standard error of a (weighted) sample quantile: sqrt(p(1-p)/ess) / f(q)
ref_quantile_se <- function(x, w, probs, q) {
  d <- stats::density(x, weights = w / sum(w), n = 2048)
  f <- stats::approx(d$x, d$y, xout = q)$y
  sqrt(probs * (1 - probs) / kish_ess(w)) / f
}

# stan_draws: matrix [iteration, chain]
# ref: either list(q = <exact quantiles>) or list(x = <draws>, w = <weights>)
compare_quantiles <- function(stan_draws, ref, label) {
  # stuck chains can contain non-finite draws; these iterations are dropped
  # here and counted in diagnostics()
  stan_draws <- stan_draws[
    apply(is.finite(stan_draws), 1, all), ,
    drop = FALSE
  ]
  q_stan <- unname(stats::quantile(stan_draws, probs))
  mcse <- unname(posterior::mcse_quantile(stan_draws, probs))
  if (!is.null(ref$q)) {
    q_ref <- ref$q
    se_ref <- rep(0, length(probs))
  } else {
    w <- if (is.null(ref$w)) rep(1, length(ref$x)) else ref$w
    q_ref <- wquantile(ref$x, w, probs)
    se_ref <- ref_quantile_se(ref$x, w, probs, q_ref)
  }
  data.frame(
    label = label,
    prob = probs,
    stan = q_stan,
    ref = q_ref,
    mcse_stan = mcse,
    se_ref = se_ref,
    z = (q_stan - q_ref) / sqrt(mcse^2 + se_ref^2)
  )
}

# prior-only fit with captured console output and warnings
fit_prior <- function(priors, cov_prior, K = 3, seed = 1,
                      iter_sampling = 5000, iter_warmup = 1000) {
  warn <- character()
  output <- utils::capture.output(
    fit <- withCallingHandlers(
      stan_gvar(ts_data[1:50, 1:K],
        priors = priors,
        cov_prior = cov_prior,
        n_chains = 4,
        iter_sampling = iter_sampling,
        iter_warmup = iter_warmup,
        refresh = 0,
        seed = seed,
        prior_only = TRUE
      ),
      warning = function(w) {
        warn <<- c(warn, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    )
  )
  list(fit = fit, output = output, warnings = warn)
}

stan_array <- function(res, par) {
  rstan::extract(res$fit$stan_fit, pars = par, permuted = FALSE)
}

diagnostics <- function(res, label) {
  sf <- res$fit$stan_fit
  summ <- suppressWarnings(rstan::summary(sf, pars = c("Beta", "Rho"))$summary)
  # diagonal of Rho is constant zero
  summ <- summ[is.finite(summ[, "Rhat"]), , drop = FALSE]
  rho <- rstan::extract(sf, "Rho")$Rho
  message("done: ", label)
  data.frame(
    label = label,
    draws = length(rstan::extract(sf, "lp__")$lp__),
    nonfinite_draws = sum(!apply(is.finite(rho), 1, all)),
    divergent = rstan::get_num_divergent(sf),
    max_treedepth = rstan::get_num_max_treedepth(sf),
    rejection_msgs = sum(grepl("reject|inverse_spd", res$output)),
    max_rhat = if (nrow(summ) > 0) max(summ[, "Rhat"]) else NA,
    min_n_eff = if (nrow(summ) > 0) min(summ[, "n_eff"]) else NA,
    n_warnings = length(res$warnings)
  )
}

pairs_upper <- function(K) {
  idx <- which(upper.tri(diag(K)), arr.ind = TRUE)
  idx[order(idx[, 1], idx[, 2]), , drop = FALSE]
}
