#' Compare two BGGM gVAR models
#'
#' @description
#' Computes the empirical distance between two models based on their point estimates
#' and compares them using reference distributions created from their
#' posterior distributions. Returns the p-value for the comparison
#' based on a decision rule specified by the user. Details are availabel in
#' TODO ADD REFERENCE TO PREPRINT.
#' @param fit_a
#' Fitted model object for Model A.
#' @param fit_b
#' Fitted model object for Model B.
#' @param cutoff
#' The percentage level of the test (default: 5\%) as integer.
#' @param dec_rule
#' The decision rule to be used. Currently only supports default "OR".
#' @param n_draws
#' The number of draws to use for reference distributions (default: 1000).
#' @param comp
#' The distance metric to use. Should be one of "frob" (Frobenius norm), "maxdiff" (maximum  difference), or "l1" (L1 norm) (default: "frob").
#' The use of the Frobenius norm is recommended.
#' @param return_all
#' Logical indicating whether to return all distributions (default: FALSE).
#' Has to be set to TRUE for plotting the results.
#' @return A list containing the results of the comparison. The list includes:
#'  \itemize{
#'   \item{sig_beta}{Binary decision on whether there is a significant difference between the temporal networks of A and B}
#'   \item{sig_pcor}{Binary decision on whether there is a significant difference between the contemporaneous networks of A and B}
#'   \item{res_beta}{The null distribution for the temporal networks for both models}
#'   \item{res_pcor}{The null distribution for the contemporaneous networks for both models}
#'   \item{emp_beta}{The empirical distance between the two temporal networks}
#'   \item{emp_pcor}{The empirical distance between the two contemporaneous networks}
#'   \item{larger_beta}{The number of reference distances larger than the empirical distance for the temporal network}
#'   \item{larger_pcor}{The number of reference distances larger than the empirical distance for the temporal network}
#'    }
#' @importFrom dplyr group_by summarize pull
#' @importFrom ggplot2 geom_density theme_classic scale_y_continuous
#' @importFrom cowplot get_legend plot_grid
#' @importFrom ggokabeito scale_fill_okabe_ito
#' @export

compare_gvar <- function(fit_a,
                         fit_b,
                         cutoff = 5,
                         dec_rule = "OR",
                         n_draws = 1000,
                         comp = "frob",
                         return_all = FALSE) {
  # Store arguments
  args <- list(match.call)

  # Input checks.
  if (!inherits(fit_a, "var_estimate")) {
    stop("Please provide a var_estimate object as input for fit_a.")
  }
  if (!inherits(fit_b, "var_estimate")) {
    stop("Please provide a var_estimate object as input for fit_b.")
  }
  # Check cutoff
  if (cutoff < 0 || cutoff > 100) {
    stop("Error: 'cutoff' must be between 0 and 100.")
  }
  if (cutoff < 1) {
    warning("Error: 'cutoff' is less than 1, which means that the alpha level is < 0.01.
            Make sure you specified the correct input.")
  }
  # Check dec_rule
  valid_dec_rules <- c("OR")
  if (!dec_rule %in% valid_dec_rules) {
    stop("Error: 'dec_rule' can only be 'OR'.")
  }

  # Check comp
  valid_comps <- c("frob", "l1", "maxdiff")
  if (!comp %in% valid_comps) {
    stop("Error: 'comp' can only be 'frob', 'l1', or 'maxdiff'.")
  }

  # Check n_draws
  if (n_draws < 1000) {
    warning("Warning: 'n_draws' below 1000 has not been tested yet.")
  }



  ## Helper function for computing distance metrics
  compute_metric <- function(a, b, metric) {
    tryCatch(
      {
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

  ## Create reference distributions for both models
  ref_a <- post_distance_within(fit_a,
    comp = comp,
    pred = FALSE,
    draws = n_draws
  )
  ref_b <- post_distance_within(fit_b,
    comp = comp,
    pred = FALSE,
    draws = n_draws
  )

  ## Empirical distance
  # Compute empirical distance as test statistic
  emp_beta <- compute_metric(fit_a$beta_mu, fit_b$beta_mu, comp)
  emp_pcor <- compute_metric(fit_a$pcor_mu, fit_b$pcor_mu, comp)


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

  ## Implement decision rule "OR"
  # Helper function
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

  if (dec_rule == "OR") {
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
      args = args
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
      args = args
    )
  }


  class(l_res) <- "compare_gvar"
  return(l_res)
}

# Plotting method
plot.compare_gvar <- function(x,
                              name_a = NULL,   # set name for a manually
                              name_b = NULL,   # set name for b manually
                              ...){

  # Input check
  if(is.null(x$null)){
    stop("Reference distributions of compare_gvar must be saved using the
         argument 'return_all'=TRUE ")
  }

  # sysfonts::font_add_google("News Cycle", "news")
  # # use showtext
  # showtext::showtext_auto()

  # Exchange names
  if(!is.null(name_a)){
    name_a <- as.character(name_a)
    x$res_beta$mod <- gsub("mod_a", name_a, x$res_beta$mod)
    x$res_pcor$mod <- gsub("mod_a", name_a, x$res_pcor$mod)
  }
  if(!is.null(name_b)){
    name_b <- as.character(name_b)
    x$res_beta$mod <- gsub("mod_b", name_b, x$res_beta$mod)
    x$res_pcor$mod <- gsub("mod_b", name_b, x$res_pcor$mod)
  }



  # Plotting
  plt_beta <- ggplot(x$res_beta,
                     aes(x = .data$null,
                         fill = .data$mod))+
    geom_density(alpha = .7)+
    theme_classic()+
    ggokabeito::scale_fill_okabe_ito()+
    geom_vline(aes(xintercept = x$emp_beta),
               col = "red", lty = 1, linewidth = .75)+
    scale_y_continuous(expand = c(0,0))+
    labs(title = "Temporal",
         y = "",
         x = "Norm Value")+
    # theme_compare()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "right")


  plt_pcor <- ggplot(x$res_pcor,
                     aes(x = .data$null,
                         fill = .data$mod))+
    geom_density(alpha = .7)+
    theme_classic()+
    ggokabeito::scale_fill_okabe_ito()+
    geom_vline(aes(xintercept = x$emp_pcor),
               col = "red", lty = 1, linewidth = .75)+
    scale_y_continuous(expand = c(0,0))+
    labs(title = "Contemporaneous",
         y = "",
         x = "Norm Value")+
    # theme_compare()+
    theme(axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = "right")

  leg <- cowplot::get_legend(plt_beta)

  # Plot
  plt_tmp <- cowplot::plot_grid(plt_beta + theme(legend.position = "none"),
                                plt_pcor + theme(legend.position = "none"))

  # Add legend
  plt <- cowplot::plot_grid(plt_tmp, leg, rel_widths = c(3, .4))
  plt

}


