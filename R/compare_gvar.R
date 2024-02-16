#' Compare two Bayesian gVAR models
#'
#' @description This function compares two Bayesian Graphical Vector
#'   Autoregressive models using matrix norms to test if the observed
#'   differences between two models is reliable. It computes the empirical
#'   distance between two models based on their point estimates and compares
#'   them using reference distributions created from their posterior
#'   distributions. Returns the p-value for the comparison based on a decision
#'   rule specified by the user. Details are available in Siepe, Kloft & Heck (2023)
#'   <doi:10.31234/osf.io/uwfjc>.
#'
#' @param fit_a Fitted model object for Model A. This can be a tsnet_fit object
#'   (obtained from [stan_gvar()]), a BGGM object (obtained from
#'   [BGGM::var_estimate()]), or extracted posterior samples (obtained from
#'   [stan_fit_convert()).
#' @param fit_b Fitted model object for Model B. This can be a tsnet_fit object
#'   (obtained from [stan_gvar()]), a BGGM object (obtained from
#'   [BGGM::var_estimate()]), or extracted posterior samples (obtained from
#'   [stan_fit_convert()).
#' @param cutoff The percentage level of the test (default: 5\%) as integer.
#' @param dec_rule The decision rule to be used. Currently supports default "or"
#'   (comparing against two reference distributions) and "comb" (combining the
#'   reference distributions). The use of "or" is recommended, as "comb" is less
#'   stable.
#' @param n_draws The number of draws to use for reference distributions
#'   (default: 1000).
#' @param comp The distance metric to use. Should be one of "frob" (Frobenius
#'   norm), "maxdiff" (maximum  difference), or "l1" (L1 norm) (default:
#'   "frob"). The use of the Frobenius norm is recommended.
#' @param return_all Logical indicating whether to return all distributions
#'   (default: FALSE). Has to be set to TRUE for plotting the results.
#' @param sampling_method Draw sequential pairs of samples from the posterior,
#'   with certain distance between them ("sequential") or randomly from two
#'   halves of the posterior ("random"). The "random" method is preferred to
#'   account for potential autocorrelation between subsequent samples. Default:
#'   "random".
#' @param indices A list of "beta" and "pcor" indices specifying which elements
#'   of the matrices to consider when calculating distances. If NULL (default),
#'   all elements of both matrices are considered. If provided, only the
#'   elements at these indices are considered. If only one of the matrices
#'   should have indices, the other one should be NULL. This can be useful if
#'   you want to calculate distances based on a subset of the elements in the
#'   matrices.
#' @param burnin The number of burn-in iterations to discard (default: 0).
#' @return A list (of class "compare_gvar") containing the results of the
#'   comparison. The list includes:
#'  \itemize{
#'   \item{sig_beta}{Binary decision on whether there is a significant difference between the temporal networks of A and B}
#'   \item{sig_pcor}{Binary decision on whether there is a significant difference between the contemporaneous networks of A and B}
#'   \item{res_beta}{The null distribution for the temporal networks for both models}
#'   \item{res_pcor}{The null distribution for the contemporaneous networks for both models}
#'   \item{emp_beta}{The empirical distance between the two temporal networks}
#'   \item{emp_pcor}{The empirical distance between the two contemporaneous networks}
#'   \item{larger_beta}{The number of reference distances larger than the empirical distance for the temporal network}
#'   \item{larger_pcor}{The number of reference distances larger than the empirical distance for the temporal network}
#'   \item{arguments}{The arguments used in the function call}
#'    }
#' @importFrom dplyr group_by summarize pull
#' @importFrom rlang .data
#'
#' @examples
#' # use internal fit data of two individuals
#' data(fit_data)
#' test_res <- compare_gvar(fit_data[[1]],
#' fit_data[[2]],
#' n_draws = 100,
#' return_all = TRUE)
#' print(test_res)
#' @export

compare_gvar <- function(fit_a,
                         fit_b,
                         cutoff = 5,
                         dec_rule = "or",
                         n_draws = 1000,
                         comp = "frob",
                         return_all = FALSE,
                         sampling_method = "random",
                         indices = NULL,
                         burnin = 0) {


  # Store arguments
  mc <- match.call()

  # Get formal arguments and their defaults
  fl <- formals()

  # Identify the missing arguments
  missing_args <- setdiff(names(fl), names(mc))

  # update with missing arguments
  missing_args <- Map(function(arg, default)
    if (!is.null(default)) mc[[arg]] <- default, missing_args, fl[missing_args])

  all_args <- c(as.list(mc), missing_args)

  # Check cutoff
  if (cutoff < 0 || cutoff > 100) {
    stop("Error: 'cutoff' must be between 0 and 100.")
  }
  if (cutoff < 1) {
    warning("Warning: 'cutoff' is less than 1, which means that the alpha level is < 0.01.
            Make sure you specified the correct input.")
  }
  # Check dec_rule
  valid_dec_rules <- c("or", "comb")
  if (!dec_rule %in% valid_dec_rules) {
    stop("Error: 'dec_rule' can only be 'or' or 'comb'.")
  }

  # Check comp
  valid_comps <- c("frob", "l1", "maxdiff")
  if (!comp %in% valid_comps) {
    stop("Error: 'comp' can only be 'frob', 'l1', or 'maxdiff'.")
  }

  # Check fit input
  # fit_a and fit_b need to either be "var_estimate" or "stanfit"
  if(!(inherits(fit_a, "var_estimate") ||
       inherits(fit_a, "tsnet_fit") ||
       inherits(fit_a, "tsnet_samples"))) {
    stop("Error: 'fit_a' must be either a 'var_estimate', 'tsnet_fit', or 'tsnet_samples' object.")
  }
  if(!(inherits(fit_b, "var_estimate") ||
       inherits(fit_b, "tsnet_fit") ||
       inherits(fit_b, "tsnet_samples"))) {
    stop("Error: 'fit_b' must be either a 'var_estimate', 'tsnet_fit', or 'tsnet_samples' object.")
  }

  # Check indices
  if (!is.null(indices)) {
    if (!is.list(indices)) {
      stop("Error: 'indices' must be a list.")
    }
    if (!all(c("beta", "pcor") %in% names(indices))) {
      stop("Error: 'indices' must contain 'beta' and 'pcor'.")
    }
    if (!is.null(indices$beta)) {
      if (!is.numeric(indices$beta)) {
        stop("Error: 'indices$beta' must be numeric.")
      }
    }
    if (!is.null(indices$pcor)) {
      if (!is.numeric(indices$pcor)) {
        stop("Error: 'indices$pcor' must be numeric.")
      }
    }
  }



  ## Input conversion
  if(inherits(fit_a, "tsnet_fit")) {
    fit_a <- stan_fit_convert(fit_a,
                                     return_params = c("beta", "pcor"))
    # randomly change array index to remove ordering of chains
    ind_a <- sample(dim(fit_a$fit$beta)[3], replace = FALSE)
    fit_a$fit$beta <- fit_a$fit$beta[,,ind_a]
    fit_a$fit$pcor <- fit_a$fit$pcor[,,ind_a]

  }
  if(inherits(fit_b, "tsnet_fit")) {
    fit_b <- stan_fit_convert(fit_b,
                                     return_params = c("beta", "pcor"))
    # randomly change array index to remove ordering of chains
    ind_b <- sample(dim(fit_b$fit$beta)[3], replace = FALSE)
    fit_b$fit$beta <- fit_b$fit$beta[,,ind_b]
    fit_b$fit$pcor <- fit_b$fit$pcor[,,ind_b]
  }




  ## Helper function for computing distance metrics
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

  ## Helper function to only use upper triangle of matrix
  ut <- function(x) {
    matrix(x[upper.tri(x, diag = FALSE)])
  }


  ## Create reference distributions for both models
  ref_a <- post_distance_within(fit_a,
    comp = comp,
    pred = FALSE,
    n_draws = n_draws,
    sampling_method = sampling_method,
    indices = indices,
    burnin = burnin
  )
  ref_b <- post_distance_within(fit_b,
    comp = comp,
    pred = FALSE,
    n_draws = n_draws,
    sampling_method = sampling_method,
    indices = indices,
    burnin = burnin
  )

  ## Empirical distance
  # Compute empirical distance as test statistic
  # provide indices if specified
  if (is.null(indices)) {
    emp_beta <- compute_metric(fit_a$beta_mu,
                               fit_b$beta_mu,
                               comp,
                               indices)
  } else {
    emp_beta <- compute_metric(fit_a$beta_mu,
                               fit_b$beta_mu,
                               comp,
                               indices$beta)
  }

  if (is.null(indices)) {
    emp_pcor <- compute_metric(ut(fit_a$pcor_mu),
                               ut(fit_b$pcor_mu),
                               comp,
                               indices)
  } else {
    emp_pcor <- compute_metric(fit_a$pcor_mu,
                               fit_b$pcor_mu,
                               comp,
                               indices$pcor)
  }


  ## Combine results
  res_beta <- data.frame(
    null = c(
      unlist(ref_a[["beta"]]),
      unlist(ref_b[["beta"]])
    ),
    mod = c(
      rep("mod_a", n_draws),
      rep("mod_b", n_draws)
    ),
    emp = rep(emp_beta, n_draws * 2),
    comp = rep(comp, n_draws * 2)
  )


  res_pcor <- data.frame(
    null = c(
      unlist(ref_a[["pcor"]]),
      unlist(ref_b[["pcor"]])
    ),
    mod = c(
      rep("mod_a", n_draws),
      rep("mod_b", n_draws)
    ),
    emp = rep(emp_pcor, n_draws * 2),
    comp = rep(comp, n_draws * 2)
  )

  ## Implement decision rule "or"
  # Helper function
  if (dec_rule == "or") {
    compute_stats <- function(data, var, cutoff, n_draws) {
      sig_decision <- data |>
        dplyr::group_by(.data$mod) |>
        dplyr::summarize(
          sum_larger =
            sum(.data$null > .data$emp)
        ) |>
        dplyr::summarize(
          sig_decision =
            sum(.data$sum_larger < cutoff * (n_draws / 100))
        ) |>
        dplyr::pull(.data$sig_decision)

      sum_larger <- data |>
        dplyr::group_by(.data$mod) |>
        dplyr::summarize(
          sum_larger =
            sum(.data$null > .data$emp)
        )

      return(list(sig_decision = sig_decision, sum_larger = sum_larger))
    }
    sig_beta <- compute_stats(res_beta, "null", cutoff, n_draws)$sig_decision
    larger_beta <- compute_stats(res_beta, "null", cutoff, n_draws)$sum_larger
    sig_pcor <- compute_stats(res_pcor, "null", cutoff, n_draws)$sig_decision
    larger_pcor <- compute_stats(res_pcor, "null", cutoff, n_draws)$sum_larger
  }
  ## Decision rule "comb"
  else if(dec_rule == "comb") {
    compute_stats <- function(data, var, cutoff, n_draws) {
      sig_decision <- data |>
        # dplyr::group_by(.data$mod) |>
        dplyr::summarize(
          sum_larger =
            sum(.data$null > .data$emp)
        ) |>
        dplyr::summarize(
          sig_decision =
            sum(.data$sum_larger < cutoff * (n_draws * 2 / 100))
        ) |>
        dplyr::pull(.data$sig_decision)

      sum_larger <- data |>
        # dplyr::group_by(.data$mod) |>
        dplyr::summarize(
          sum_larger =
            sum(.data$null > .data$emp)
        )

      return(list(sig_decision = sig_decision, sum_larger = sum_larger))
    }
    sig_beta <- compute_stats(res_beta, "null", cutoff, n_draws)$sig_decision
    larger_beta <- compute_stats(res_beta, "null", cutoff, n_draws)$sum_larger
    sig_pcor <- compute_stats(res_pcor, "null", cutoff, n_draws)$sig_decision
    larger_pcor <- compute_stats(res_pcor, "null", cutoff, n_draws)$sum_larger

  }

  if (!return_all) {
    l_res <- list(
      sig_beta = sig_beta,
      sig_pcor = sig_pcor,
      emp_beta = emp_beta,
      emp_pcor = emp_pcor,
      larger_beta = larger_beta,
      larger_pcor = larger_pcor,
      arguments = all_args
    )
  }
  if (isTRUE(return_all)) {
    l_res <- list(
      sig_beta = sig_beta,
      sig_pcor = sig_pcor,
      res_beta = res_beta,
      res_pcor = res_pcor,
      emp_beta = emp_beta,
      emp_pcor = emp_pcor,
      larger_beta = larger_beta,
      larger_pcor = larger_pcor,
      arguments = all_args
    )
  }


  class(l_res) <- c("compare_gvar", class(l_res))
  return(l_res)
}


