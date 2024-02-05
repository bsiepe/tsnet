#' Plot compare_gvar
#'
#' This function is a plotting method for the class produced by
#' [compare_gvar()]. It generates a plot showing the density of posterior
#' uncertainty distributions for distances and the empirical distance value for two GVAR models.
#'
#' @param x An object of class "compare_gvar".
#' @param name_a Optional. The name for model A. If provided, it replaces
#'   "mod_a" in the plot.
#' @param name_b Optional. The name for model B. If provided, it replaces
#'   "mod_b" in the plot.
#' @param ... Additional arguments to be passed to the plotting functions.
#'
#' @details The function first checks if the full reference distributions of
#' [compare_gvar()] are saved using the argument 'return_all' set to TRUE. If
#' not, an error is thrown.
#'
#' Using the "name_a" and "name_b" arguments allows for custom labeling of the
#' two models in the plot.
#'
#' The function generates two density plots using `ggplot2`, one for the
#' temporal network (beta) and another for the contemporaneous network (pcor).
#' The density distributions are filled with different colors based on the
#' corresponding models (mod_a and mod_b). The empirical distances between the
#' networks are indicated by red vertical lines.
#'
#' @import ggplot2
#' @importFrom ggokabeito scale_fill_okabe_ito palette_okabe_ito
#' @importFrom cowplot plot_grid get_legend
#' @importFrom ggdist theme_ggdist
#'
#' @examples 
#' \dontrun{
#' data(fit_data)
#' test_res <- compare_gvar(fit_data[[1]], fit_data[[2]])
#' plot(test_res)
#' }
#' @export
plot.compare_gvar <- function(x,
                              name_a = NULL,
                              name_b = NULL,
                              ...){

  # Input check
  if(is.null(x$res_beta$null) & is.null(x$res_pcor$null)){
    stop("Reference distributions of compare_gvar must be saved using the
         argument 'return_all'=TRUE ")
  }
  if(!inherits(x, "compare_gvar")){
    stop("This function only works with a result of the compare_gvar function.")
  }


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

  # Check test type
  if(x$arguments$dec_rule == "or"){
    # Plotting
    plt_beta <- ggplot2::ggplot(x$res_beta,
                                aes(x = .data$null,
                                    fill = .data$mod))+
      ggplot2::geom_density(alpha = .7)+
      ggdist::theme_ggdist() +
      ggokabeito::scale_fill_okabe_ito(order = c(5, 1:8))+
      ggplot2::geom_vline(aes(xintercept = x$emp_beta),
                          col = "black", lty = 1, linewidth = .75)+
      ggplot2::scale_y_continuous(expand = c(0,0))+
      ggplot2::labs(title = "Temporal",
                    y = "",
                    x = "Norm Value",
                    fill = "Model")+

      ggplot2::theme(axis.ticks.y = element_blank(),
                     axis.text.y = element_blank(),
                     legend.position = "right")


    plt_pcor <- ggplot2::ggplot(x$res_pcor,
                                aes(x = .data$null,
                                    fill = .data$mod))+
      ggplot2::geom_density(alpha = .7)+
      ggdist::theme_ggdist() +
      ggokabeito::scale_fill_okabe_ito(order = c(5, 1:8))+
      ggplot2::geom_vline(aes(xintercept = x$emp_pcor),
                          col = "black", lty = 1, linewidth = .75)+
      ggplot2::scale_y_continuous(expand = c(0,0))+
      ggplot2::labs(title = "Contemporaneous",
                    y = "",
                    x = "Norm Value",
                    fill = "Model")+
      ggplot2::theme(axis.ticks.y = element_blank(),
                     axis.text.y = element_blank(),
                     legend.position = "right")

    leg <- cowplot::get_legend(plt_beta)

    # Plot
    plt_tmp <- cowplot::plot_grid(plt_beta + theme(legend.position = "none"),
                                  plt_pcor + theme(legend.position = "none"))

    # Add legend
    plt <- cowplot::plot_grid(plt_tmp, leg, rel_widths = c(3, .4))
  }
  else if(x$arguments$dec_rule == "comb"){
    # Plotting
    plt_beta <- ggplot2::ggplot(x$res_beta,
                                aes(x = .data$null))+
      ggplot2::geom_density(alpha = .7,
                            fill = ggokabeito::palette_okabe_ito()[5])+
      ggdist::theme_ggdist() +
      ggplot2::geom_vline(aes(xintercept = x$emp_beta),
                          col = "black", lty = 1, linewidth = .75)+
      ggplot2::scale_y_continuous(expand = c(0,0))+
      ggplot2::labs(title = "Temporal",
                    y = "",
                    x = "Norm Value",
                    fill = "Model")+

      ggplot2::theme(axis.ticks.y = element_blank(),
                     axis.text.y = element_blank(),
                     legend.position = "right")


    plt_pcor <- ggplot2::ggplot(x$res_pcor,
                                aes(x = .data$null))+
      ggplot2::geom_density(alpha = .7,
                            fill = ggokabeito::palette_okabe_ito()[5])+
      ggdist::theme_ggdist() +
      ggplot2::geom_vline(aes(xintercept = x$emp_pcor),
                          col = "black", lty = 1, linewidth = .75)+
      ggplot2::scale_y_continuous(expand = c(0,0))+
      ggplot2::labs(title = "Contemporaneous",
                    y = "",
                    x = "Norm Value",
                    fill = "Model")+
      ggplot2::theme(axis.ticks.y = element_blank(),
                     axis.text.y = element_blank(),
                     legend.position = "right")


    # Plot
    plt <- cowplot::plot_grid(plt_beta + theme(legend.position = "none"),
                                  plt_pcor + theme(legend.position = "none"))

  }

  print(plt)
  
}

