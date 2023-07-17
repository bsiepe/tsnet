#' Print method for compare_gvar objects
#'
#' This function prints a summary of the Norm-Based Comparison Test for a [compare_gvar()] object.
#'
#' @param x A compare_gvar object.
#' @param ... Additional arguments to be passed to the print method. (currently not used)
#'
#' @return None
#'
#' @details This function prints a summary of the Norm-Based Comparison Test for a compare_gvar object.
#' It displays the general summary and model-specific results, including the number of significant comparisons
#' in the temporal and contemporaneous networks, as well as the number of reference distances that were larger
#' than the empirical distance for each network.
#'
#'
#' @export
print.compare_gvar <- function(x,
                                 ...){

  if(!inherits(x, "compare_gvar")){
    stop("This function only works with a result of the compare_gvar function.")
  }
  cat("### Summary of the Norm-Based Comparison Test ###")
  cat("\n")
  cat("\n#--- General Summary ---#")
  cat(
    "\nIn the temporal network", x$sig_beta, "of the 2 comparisons were significant."
  )
  cat(
    "\nIn the contemporaneous network", x$sig_pcor, "of the 2 comparisons were significant."
  )
  cat("\n")
  cat(
    "\n#--- Model-specific Results ---#"
  )
  cat(
    "\nFor", x$larger_beta$mod[1], x$larger_beta$sum_larger[1],
    "of the reference distances of the temporal network and",
    x$larger_pcor$sum_larger[1],
    "of the reference distances of the contemporaneous network were larger than the empirical distance."
  )
  cat("\n")
  cat(
    "\nFor", x$larger_beta$mod[2], x$larger_beta$sum_larger[2],
    "of the reference distances of the temporal network and",
    x$larger_pcor$sum_larger[2],
    "of the reference distances of the contemporaneous network were larger than the empirical distance."
  )

}
