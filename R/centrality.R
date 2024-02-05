#' Compute Centrality Measures
#'
#' This function computes various network centrality measures for a given GVAR
#' fit object. Centrality measures describe the "connectedness" of a variable in a network, while density describes the networks' overall connectedness.
#' Specifically, it computes the in-strength, out-strength,
#' contemporaneous strength, temporal network density, and contemporaneous
#' network density. The result can then be visualized using [plot_centrality()].
#'
#' @param fitobj
#' Fitted model object for a Bayesian GVAR model.
#' This can be a stanfit object (obtained from [stan_gvar()]),
#' a BGGM object (obtained from [BGGM::var_estimate()]),
#' or extracted posterior samples (obtained from [stan_fit_convert()).
#' @param burnin An integer specifying the number of initial samples to discard
#'   as burn-in. Default is 0.
#' @param remove_ar A logical value specifying whether to remove the
#'   autoregressive effects for centrality calculation. Default is TRUE. This is
#'   only relevant for the calculation of temporal centrality/density measures.
#'
#' @return A list containing the following centrality measures:
#' \itemize{
#'   \item \code{instrength}: In-strength centrality.
#'   \item \code{outstrength}: Out-strength centrality.
#'   \item \code{strength}: Contemporaneous strength centrality.
#'   \item \code{density_beta}: Temporal network density.
#'   \item \code{density_pcor}: Contemporaneous network density.
#' }
#'
#' @examples
#'  # Use first individual from example fit data from tsnet
#'  data(fit_data)
#'  centrality_measures <- get_centrality(fit_data[[1]])
#'
#' @export
get_centrality <- function(fitobj,
                           burnin = 0,
                           remove_ar = TRUE) {

  # Input checks
  if(!(inherits(fitobj, "var_estimate") ||
       inherits(fitobj, "stanfit") ||
       inherits(fitobj, "tsnet_samples"))) {
    stop("Error: 'fit_a' must be either a 'var_estimate', 'stanfit', or 'tsnet_samples' object.")
  }

  # Input Conversion
  if(inherits(fitobj, "stanfit")) {
    fit_a <- tsnet::stan_fit_convert(fit_a,
                                     return_params = c("beta", "pcor"))
  }


  # Obtain samples
  n_samps <- dim(fitobj$fit$beta)[3]

  beta_samps <- abs(fitobj$fit$beta[, , burnin:n_samps])
  pcor_samps <- abs(fitobj$fit$pcors[, , burnin:n_samps])

  cnames <- colnames(fitobj$Y)

  if(isTRUE(remove_ar)){
    # Function to set the diagonal elements of a matrix to zero
    diag_zero <- function(mat) {
      diag(mat) <- 0
      return(mat)
    }
    beta_samps <- array(apply(beta_samps, MARGIN = 3, FUN = diag_zero),
                        dim = dim(beta_samps))
  }

  #--- Centrality measures
  # In-strength
  instrength <- t(apply(beta_samps, MARGIN = 3, FUN = colSums))
  colnames(instrength) <- cnames

  # Out-strength
  outstrength <- t(apply(beta_samps, MARGIN = 3, FUN = rowSums))
  colnames(outstrength) <- cnames

  # Contemporaneous strength
  strength <- t(apply(pcor_samps, MARGIN = 3, FUN = rowSums))
  colnames(strength) <- cnames

  # Density
  density_beta <- apply(beta_samps, MARGIN = 3, FUN = sum)
  density_pcor <- apply(pcor_samps, MARGIN = 3, FUN = sum)

  #--- Return
  return(list(
    instrength = instrength,
    outstrength = outstrength,
    strength = strength,
    density_beta = density_beta,
    density_pcor = density_pcor
  ))
}



#' Plot Centrality Measures
#'
#' This function creates a plot of various centrality measures for a given
#' object. The plot can be either a "tiefighter" plot or a "density" plot. The
#' "tiefighter" plot shows the centrality measures for each variable with
#' uncertainty bands, while the "density" plot shows the full density of the
#' centrality measures.
#'
#' @param obj An object containing the centrality measures obtained from
#'   [get_centrality()].
#' @param plot_type A character string specifying the type of plot. Accepts
#'   "tiefighter" or "density". Default is "tiefighter".
#' @param cis A numeric value specifying the credible interval. Must be between
#'   0 and 1 (exclusive). Default is 0.95.
#'
#' @return A ggplot object visualizing the centrality measures. For a
#'   "tiefighter" plot, each point represents the mean centrality measure for a
#'   variable, and the bars represent the credible interval. In a "density"
#'   plot, distribution of the centrality measures is visualized.
#'
#' @examples
#' \dontrun{
#' data(fit_data)
#' obj <- get_centrality(fit_data[[1]])
#'   plot_centrality(obj,
#'   plot_type = "tiefighter",
#'   cis = 0.95)
#' }
#' 
#' @importFrom tidyr pivot_longer
#' @importFrom dplyr group_by summarize
#' @importFrom stats quantile
#' @import ggplot2
#' 
#' @export
plot_centrality <- function(obj,
                            plot_type = "tiefighter",
                            cis = 0.95) {
  if (!is.numeric(cis) || any(cis <= 0) || any(cis >= 1)) {
    stop("cis must be a numeric vector with values between 0 and 1 (exclusive)")
  }


  #--- Preparation
  # Combine all centrality measures
  create_centrality_df <- function(centrality, suffix) {
    df <- as.data.frame(centrality)
    colnames(df) <- paste0(colnames(df), "_", suffix)
    return(df)
  }

  instrength <- create_centrality_df(obj$instrength, "Temporal\nInstrength")
  outstrength <- create_centrality_df(obj$outstrength, "Temporal\nOutstrength")
  strength <- create_centrality_df(obj$strength, "Contemporaneous\nStrength")


  df_centrality <- cbind(
    instrength,
    outstrength,
    strength
  )

  #--- Overview
  if (plot_type == "tiefighter") {
    #--- Plot
    overview_plot <- df_centrality %>%
      tidyr::pivot_longer(
        cols = everything(),
        names_to = "measure",
        values_to = "value"
      ) |>
      dplyr::group_by(.data$measure) |>
      dplyr::summarize(
        mean_value = mean(.data$value,
                          na.rm = TRUE
        ),
        lb = stats::quantile(.data$value,
                      probs = (1 - cis) / 2,
                      na.rm = TRUE
        ),
        ub = stats::quantile(.data$value,
                      probs = 1 - ((1 - cis) / 2)
        ),
        na.rm = TRUE
      ) |>
      dplyr::ungroup() |>
      tidyr::separate(
        col = .data$measure,
        into = c("variable", "centrality"),
        sep = "_"
      ) |>
      ggplot2::ggplot(aes(x = .data$mean_value, y = .data$variable)) +
      ggplot2::geom_point() +
      ggplot2::geom_errorbarh(aes(
        xmin = .data$lb,
        xmax = .data$ub
      )) +
      ggdist::theme_ggdist() +
      ggplot2::facet_wrap(. ~ .data$centrality) +
      ggplot2::labs(
        x = "Centrality",
        y = "Variable"
      )

    return(overview_plot)
  }
  #--- Density
  if (plot_type == "density") {
    #--- Plot
    density_plot <- df_centrality %>%
      tidyr::pivot_longer(
        cols = everything(),
        names_to = "measure",
        values_to = "value"
      ) |>
      dplyr::group_by(.data$measure) |>
      dplyr::mutate(mean_value = mean(.data$value,
                                      na.rm = TRUE
      )) |>
      dplyr::ungroup() |>
      tidyr::separate(
        col = .data$measure,
        into = c("variable", "centrality"),
        sep = "_"
      ) |>
      ggplot2::ggplot(aes(
        x = .data$value,
        y = .data$variable
      )) +
      ggdist::stat_slab(
        aes(
          fill = after_stat(.data$level)
        ),
        # fixed width for now
        .width = c(0.8, 0.9, 0.95, 1)
      ) +
      ggdist::stat_pointinterval(aes(x = .data$mean_value),
                                 size = 1
      ) +
      ggplot2::scale_alpha(guide = "none") +
      ggdist::theme_ggdist() +
      ggplot2::facet_wrap(. ~ .data$centrality) +
      ggplot2::scale_fill_brewer() +
      ggplot2::labs(
        x = "Centrality",
        y = "Variable",
        fill = "Credible Interval"
      )

    return(density_plot)
  }
}

