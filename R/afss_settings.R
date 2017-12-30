#' Generate a list of settings for automated factor slice sampling
#'
#' @param factor_update_interval number of iterations between factor and
#'   interval updates (defaults to an update every iteration)
#' @param initial_widths vector of initial slice widths
#' @param initial_factors matrix of initial factor loadings
#'
#' @return list with additional settings for automated factor slice sampling
#' @export
afss_settings <- function(factor_update_interval = 1, initial_widths = NULL, initial_factors = NULL) {
      
      if((is.null(initial_widths) & !is.null(initial_factors)) | 
         (!is.null(initial_widths) & is.null(initial_factors))) {
            stop("If one of initial widths or factors are supplied, the other must be supplied.")
      }
      
      list(factor_update_interval = factor_update_interval,
           initial_widths = initial_widths,
           initial_factors = initial_factors)
}