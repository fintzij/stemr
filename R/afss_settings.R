#' Generate a list of settings for automated factor slice sampling
#'
#' @param factor_update_interval number of iterations between factor and
#'   interval updates (defaults to an update every iteration)
#' @param weight_update_interval number of iterations between updates to the
#'   factor sampling weights, defaults to the factor update interval.
#' @param n_fss_updates number of slice directions that are sampled per
#'   iteration, defaults to the number of model parameters if NULL.
#' @param initial_widths vector of initial slice widths
#' @param initial_factors matrix of initial factor loadings
#' @param adapt_factors should the factors be updated? defaults to TRUE
#' @param adapt_intervals should the interval 
#'
#' @return list with additional settings for automated factor slice sampling
#' @export
afss_settings <-
      function(factor_update_interval = 1,
               weight_update_interval = NULL,
               n_fss_updates = NULL,
               adapt_factors = TRUE,
               initial_widths = NULL,
               initial_factors = NULL) {
            
      if((is.null(initial_widths) & !is.null(initial_factors)) | 
         (!is.null(initial_widths) & is.null(initial_factors))) {
            stop("If one of initial widths or factors are supplied, the other must be supplied.")
      }
      
      if(is.null(weight_update_interval)) weight_update_interval <- factor_update_interval
      
      list(factor_update_interval = factor_update_interval,
           weight_update_interval = weight_update_interval,
           n_fss_updates          = n_fss_updates,
           adapt_factors          = adapt_factors,
           initial_widths         = initial_widths,
           initial_factors        = initial_factors)
}