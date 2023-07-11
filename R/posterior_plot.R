#' posterior_plot
#'
#' Plots posterior distributions of the parameters of the temporal and
#' contemporaneous networks of a gVAR model.
#'
#' @param fitobj A 'var_estimate' fit object from the `BGGM` package.
#' @param mat A matrix to use for plotting. Possibilities include "beta" (temporal network)
#' and "pcor" (contemporaneous network). Default is "beta" (temporal network).
#' @param cis A numeric vector of credible intervals to use for plotting. Default is c(0.8, 0.9, 0.95).
#'
#' @import ggdist
#' @import tidyr
#' @import dplyr
#' @importFrom BGGM posterior_samples
#' @importFrom ggplot2 aes facet_grid geom_vline labs scale_alpha scale_fill_brewer ylim theme element_blank ggplot after_stat
#' @importFrom ggdist stat_pointinterval stat_slab
#' @importFrom tidyr separate_wider_delim pivot_longer
#' @export

posterior_plot <- function(fitobj,
                           mat = "beta",
                           cis = c(0.8, 0.9, 0.95)) { # credible intervals for plotting
  # Input Checks
  if (!inherits(fitobj, "var_estimate")) {
    stop("Please provide a var_estimate object as input for fit_a.")
  }

  if (!is.numeric(cis) || any(cis <= 0) || any(cis >= 1)) {
    stop("cis must be a numeric vector with values between 0 and 1 (exclusive)")
  }

  if (length(grep("_", colnames(fitobj$Y))) > 0) {
    stop("Column names must not contain an underscore. Please rename.")
  }


  # Obtain samples
  samps <- BGGM::posterior_samples(fitobj)


  # Split into betas and pcors
  beta_cols <- grep(".l1", colnames(samps), value = TRUE)
  pcor_cols <- grep("--", colnames(samps), value = TRUE)

  beta_samps <- as.data.frame(samps[, beta_cols])
  pcor_samps <- as.data.frame(samps[, pcor_cols])

  # Pivot longer
  beta <- beta_samps |>
    as.data.frame() |>
    dplyr::mutate(iteration = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = !.data$iteration, names_to = "edge", values_to = "value") |>
    # split edge description into nodes
    tidyr::separate_wider_delim(
      cols = .data$edge, delim = "_",
      names = c("dv", "iv")
    )

  pcor <- pcor_samps |>
    as.data.frame() |>
    dplyr::mutate(iteration = dplyr::row_number()) |>
    tidyr::pivot_longer(cols = !.data$iteration, names_to = "edge", values_to = "value") |>
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
      ggplot(aes(x = .data$value)) +
      # ggdist::stat_halfeye(aes(fill = after_stat(level)), .width = cis)+
      ggdist::stat_slab(aes(fill = after_stat(.data$level), alpha = abs(.data$mean_value)), .width = c(cis, 1)) +
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
      labs(
        y = "",
        fill = "CI"
      ) +
      ylim(-0.1, 1)


    print(beta_plot)
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
      ggplot(aes(x = .data$value)) +
      ggdist::stat_slab(aes(fill = after_stat(.data$level), alpha = abs(.data$mean_value)), .width = c(cis, 1)) +
      ggdist::stat_pointinterval(aes(alpha = abs(.data$mean_value)), size = 1) +
      facet_grid(iv ~ dv,
        switch = "y"
      ) +
      ggdist::theme_ggdist() +
      scale_alpha(guide = "none") +
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

    print(pcor_plot)
  }
}
