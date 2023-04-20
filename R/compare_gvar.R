#' Compare Variance Components between Two Models
#'
#' @description
#' Computes the empirical distance between two models based on the variance components
#' and compares them using reference distributions. Returns the p-value for the comparison
#' based on a decision rule specified by the user.
#' @param fit_a
#' Fitted model object for Model A
#' @param fit_b
#' Fitted model object for Model B
#' @param cutoff
#' The percentage level of the test (default: 5\%)
#' @param dec_rule
#' The decision rule to be used. Currently only supports default "OR".
#' @param n_draws
#' The number of draws to use for reference distributions (default: 1000)
#' @param comp
#' The distance metric to use. Should be one of "frob" (Frobenius norm), "maxdiff" (maximum  difference), or "l1" (L1 norm) (default: "frob")
#' @param return_all
#' Logical indicating whether to return all distributions (default: FALSE)
#' @return A list containing the results of the comparison. The list includes:
#'  \itemize{
#'   \item{sig_beta}{The decision on whether there is a significant difference between the variance components for Model A and Model B (based on the beta parameter)}
#'   \item{sig_pcor}{The decision on whether there is a significant difference between the variance components for Model A and Model B (based on the partial correlation parameter)}
#'   \item{res_beta}{The null distribution for the variance components (based on the beta parameter) for both models}
#'   \item{res_pcor}{The null distribution for the variance components (based on the partial correlation parameter) for both models}
#'   \item{emp_beta}{The empirical distance between the two models (based on the beta parameter)}
#'   \item{emp_pcor}{The empirical distance between the two models (based on the partial correlation parameter)}
#'   \item{larger_beta}{The number of times the null hypothesis (based on the beta parameter) was rejected across all draws}
#'   \item{larger_pcor}{The number of times the null hypothesis (based on the partial correlation parameter) was rejected across all draws}
#'    }
#' @importFrom dplyr group_by summarize pull
#' @export

compare_gvar <- function(fit_a,
                        fit_b,
                        cutoff = 5,
                        dec_rule = "OR",
                        n_draws = 1000,
                        comp = "frob",
                        return_all = FALSE){


  # Store arguments
  args <- list(match.call)

  ## Helper function for computing distance metrics
  compute_metric <- function(a, b, metric) {
    tryCatch({
      if (metric == "frob") {
        norm(a - b, type = "F")
      } else if (metric == "maxdiff") {
        max(abs(a - b))
      } else if (metric == "l1") {
        sum(abs(a - b))
      }
    }, error = function(e) NA)
  }

  ## Create reference distributions for both models
  ref_a <- post_distance_within(fit_a,
                                comp = comp,
                                pred = FALSE,
                                draws = n_draws)
  ref_b <- post_distance_within(fit_b,
                                comp = comp,
                                pred = FALSE,
                                draws = n_draws)

  ## Empirical distance
  # Compute empirical distance as test statistic
  emp_beta <- compute_metric(fit_a$beta_mu, fit_b$beta_mu, comp)
  emp_pcor <- compute_metric(fit_a$pcor_mu, fit_b$pcor_mu, comp)


  ## Combine results
  res_beta <- data.frame(null = c(unlist(ref_a[["beta"]]),
                                  unlist(ref_b[["beta"]])),
                         mod = c(rep("mod_a", n_draws),
                                 rep("mod_b", n_draws)),
                         emp = rep(emp_beta, n_draws*2),
                         comp = rep(comp, n_draws*2))


  res_pcor <- data.frame(null = c(unlist(ref_a[["pcor"]]),
                                  unlist(ref_b[["pcor"]])),
                         mod = c(rep("mod_a", n_draws),
                                 rep("mod_b", n_draws)),
                         emp = rep(emp_pcor, n_draws*2),
                         comp = rep(comp, n_draws*2))

  ## Implement decision rule "OR"
  # Helper function
  compute_stats <- function(data, var, cutoff, n_draws) {
    sig_decision <- data |>
      dplyr::group_by(.data$mod) |>
      dplyr::summarize(sum_larger =
                         sum(.data$null > .data$emp)) |>
      dplyr::summarize(sig_decision =
                         sum(.data$sum_larger < cutoff * (n_draws/100))) |>
      dplyr::pull(.data$sig_decision)

    sum_larger <- data |>
      dplyr::group_by(.data$mod) |>
      dplyr::summarize(sum_larger =
                         sum(.data$null > .data$emp))

    return(list(sig_decision = sig_decision, sum_larger = sum_larger))
  }

  if(dec_rule == "OR"){
    sig_beta <- compute_stats(res_beta, "null", cutoff, n_draws)$sig_decision
    larger_beta <- compute_stats(res_beta, "null", cutoff, n_draws)$sum_larger
    sig_pcor <- compute_stats(res_pcor, "null", cutoff, n_draws)$sig_decision
    larger_pcor <- compute_stats(res_pcor, "null", cutoff, n_draws)$sum_larger

  }

  if(!return_all){
    l_res <- list(sig_beta = sig_beta,
                  sig_pcor = sig_pcor,
                  emp_beta = emp_beta,
                  emp_pcor = emp_pcor,
                  larger_beta = larger_beta,
                  larger_pcor = larger_pcor,
                  args = args)

  }
  if(isTRUE(return_all)){
    l_res <- list(sig_beta = sig_beta,
                  sig_pcor = sig_pcor,
                  res_beta = res_beta,
                  res_pcor = res_pcor,
                  emp_beta = emp_beta,
                  emp_pcor = emp_pcor,
                  larger_beta = larger_beta,
                  larger_pcor = larger_pcor,
                  args = args)

  }



  return(l_res)



}
