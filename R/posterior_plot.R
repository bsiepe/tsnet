#' posterior_plot
#'
#' @description Plots posterior distributions of the parameters of the temporal
#' or the contemporaneous networks of a GVAR model. The posterior distributions
#' are visualized as densities in a matrix layout.
#'
#' @param fitobj Fitted model object. This can be a tsnet_fit object (obtained
#'   from [stan_gvar()]) or a BGGM object (obtained from [BGGM::var_estimate()]).
#' @param mat A matrix to use for plotting. Possibilities include "beta"
#'   (temporal network) and "pcor" (contemporaneous network). Default is "beta"
#'   (temporal network).
#' @param cis A numeric vector of credible intervals to use for plotting.
#'   Default is c(0.8, 0.9, 0.95).
#'
#' @details In the returned plot, posterior distributions for every parameter
#' are shown. Lagged variables are displayed along the vertical line of the
#' grid, and non-lagged variables along the horizontal line of the grids.
#'
#' @import ggplot2
#' @importFrom ggdist stat_pointinterval stat_slab
#' @importFrom tidyr separate_wider_delim pivot_longer
#'
#' @examples
#' \dontrun{
#' # Load simulated time series data
#' data(ts_data)
#' example_data <- ts_data[1:100,1:4]
#'
#' # Estimate a GVAR model
#' fit <- stan_gvar(example_data, n_chains = 2)
#'
#' # Extract posterior samples
#' posterior_plot(fit)
#' }
#' @export

posterior_plot <- function(fitobj,
                           mat = "beta",
                           cis = c(0.8, 0.9, 0.95)) {

  # browser()

  if (!is.numeric(cis) || any(cis <= 0) || any(cis >= 1)) {
    stop("cis must be a numeric vector with values between 0 and 1 (exclusive)")
  }


  samps <- prepare_samples_plot(fitobj)




  # Split into betas and pcors
  beta_cols <- grep(".l1", colnames(samps), value = TRUE)
  pcor_cols <- grep("--", colnames(samps), value = TRUE)

  beta_samps <- as.data.frame(samps[, beta_cols])
  pcor_samps <- as.data.frame(samps[, pcor_cols])

  # order of variables for plotting
  if(inherits(fitobj, "tsnet_fit")) {
    if (length(grep("_", colnames(fitobj$arguments$cnames))) > 0) {
      stop("Column names must not contain an underscore. Please rename.")
    }
    beta_order <- paste0(fitobj$arguments$cnames, ".l1")
    pcor_order <- fitobj$arguments$cnames
  } else {
    if (length(grep("_", colnames(fitobj$Y))) > 0) {
      stop("Column names must not contain an underscore. Please rename.")
    }
    beta_order <- rownames(fitobj$beta_mu)
    pcor_order <- colnames(fitobj$beta_mu)
  }



  # Pivot longer
  beta <- beta_samps |>
    as.data.frame() |>
    dplyr::mutate(iteration = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = !.data$iteration,
                        names_to = "edge",
                        values_to = "value") |>
    # split edge description into nodes
    tidyr::separate_wider_delim(
      cols = .data$edge, delim = "_",
      names = c("dv", "iv")
    )

  pcor <- pcor_samps |>
    as.data.frame() |>
    dplyr::mutate(iteration = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = !.data$iteration,
                        names_to = "edge",
                        values_to = "value") |>
    # split edge description into nodes
    tidyr::separate_wider_delim(
      cols = .data$edge, delim = "--",
      names = c("dv", "iv")
    )



  # Create matrix layout
  if (mat == "beta") {
    # Start plotting
    beta_plot <- beta |>
      dplyr::group_by(.data$dv, .data$iv) |>
      dplyr::mutate(mean_value = mean(.data$value, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      dplyr::mutate(dv = factor(.data$dv, levels = pcor_order),
                    iv = factor(.data$iv, levels = beta_order)) |>
      ggplot(aes(x = .data$value)) +
      # ggdist::stat_halfeye(aes(fill = after_stat(level)), .width = cis)+
      ggdist::stat_slab(aes(fill = after_stat(.data$level),
                            alpha = abs(.data$mean_value)),
                        .width = c(cis, 1)) +
      ggdist::stat_pointinterval(aes(alpha = abs(.data$mean_value)), size = 1) +
      scale_alpha(guide = "none") +
      facet_grid(iv ~ dv,
        switch = "y"
      ) +
      ggdist::theme_ggdist() +
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      scale_fill_brewer() +
      scale_x_continuous(breaks = c(-1, 0, 1),
                         minor_breaks = c(-0.5, 0.5))+
      labs(
        y = "",
        fill = "CI"
      ) +
      ylim(-0.1, 1)

      return(beta_plot)

  }

  if (mat == "pcor") {
    # make symmetric by splitting
    pcor_tmp1 <- pcor |>
      dplyr::group_by(.data$dv, .data$iv) |>
      dplyr::mutate(mean_value = mean(.data$value, na.rm = TRUE)) |>
      dplyr::ungroup()
    pcor_tmp2 <- pcor |>
      dplyr::group_by(.data$dv, .data$iv) |>
      dplyr::mutate(mean_value = mean(.data$value, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      dplyr::mutate(dv2 = .data$dv, iv2 = .data$iv) |>
      dplyr::mutate(dv = .data$iv2, iv = .data$dv2) |>
      dplyr::select(-c(.data$iv2, .data$dv2))
    pcor <- rbind(pcor_tmp1, pcor_tmp2)

    pcor_plot <- pcor |>
      dplyr::group_by(.data$dv, .data$iv) |>
      dplyr::mutate(mean_value = mean(.data$value, na.rm = TRUE)) |>
      dplyr::ungroup() |>
      dplyr::mutate(dv = factor(.data$dv, levels = pcor_order),
                    iv = factor(.data$iv, levels = pcor_order)) |>
      ggplot(aes(x = .data$value)) +
      ggdist::stat_slab(aes(fill = after_stat(.data$level),
                            alpha = abs(.data$mean_value)),
                        .width = c(cis, 1)) +
      ggdist::stat_pointinterval(aes(alpha = abs(.data$mean_value)),
                                 size = 1) +
      facet_grid(iv ~ dv,
        switch = "y"
      ) +
      ggdist::theme_ggdist() +
      scale_alpha(guide = "none") +
      scale_x_continuous(breaks = c(-1, 0, 1),
                         minor_breaks = c(-0.5, 0.5))+
      geom_vline(xintercept = 0, linetype = "dashed") +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      ) +
      scale_fill_brewer() +
      labs(
        y = "",
        fill = "CI",
        alpha = ""
      ) +
      ylim(-0.1, 1)

      return(pcor_plot)

  }
}
