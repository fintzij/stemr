#' Generate a list of settings for automated factor slice sampling
#'
#' @param factor_update_interval number of iterations between factor and
#'   interval updates (defaults to an update every iteration). May also be
#'   supplied as a function that takes the current iteration as an input and
#'   returns the next iteration at which the factors should be updated. Must
#'   have a default set so that the first update can be determined. e.g.,
#'   \code{function(x = 10) 2*x} will first update the factors at iteration 10
#'   and will subsequently update the factors at iterations 20, 40, 80, ...
#' @param weight_update_interval number of iterations between updates to the
#'   factor sampling weights, defaults to the factor update interval. May also
#'   be supplied as a function that takes the current iteration as an input and
#'   returns the next iteration at which the factors should be updated. Must
#'   have a default set so that the first update can be determined.
#' @param n_fss_updates number of slice directions that are sampled per
#'   iteration, defaults to the number of model parameters if NULL.
#' @param initial_widths vector of initial slice widths
#' @param initial_factors matrix of initial factor loadings
#' @param first_factor_update iteration at which the first factor update should
#'   occur
#' @param first_weight_update iteration at which the first update to the slice
#'   weights should occur
#' @param sample_all_initially should all factors be sampled per iteration until
#'   the first factor update, defaults to TRUE
#' @param adapt_factors should the factors be updated? defaults to TRUE
#'
#' @return list with additional settings for automated factor slice sampling
#' @export
afss_settings <-
      function(factor_update_interval = 1,
               weight_update_interval = NULL,
               first_factor_update = 100,
               first_weight_update = 100,
               sample_all_initially = TRUE,
               n_fss_updates = NULL,
               adapt_factors = TRUE,
               initial_widths = NULL,
               initial_factors = NULL) {
            
      if((is.null(initial_widths) & !is.null(initial_factors)) | 
         (!is.null(initial_widths) & is.null(initial_factors))) {
            stop("If one of initial widths or factors are supplied, the other must be supplied.")
      } 
         
      # test if the factor and weight update interval functions have defaults if they are specified via functions   
      if(is.function(factor_update_interval)) {
            fui <- NULL
            try({fui <- factor_update_interval()}, silent = T)
            if(is.null(fui)) stop("If the factor update interval is supplied as a function, it must have a default value for the first interval.")
      }
      
      if(is.function(weight_update_interval)) {
            fui <- NULL
            try({fui <- weight_update_interval()}, silent = T)
            if(is.null(fui)) stop("If the weight update interval is supplied as a function, it must have a default value for the first interval.")
      }
      
      if(is.null(weight_update_interval)) weight_update_interval <- factor_update_interval
      
      list(factor_update_interval = factor_update_interval,
           weight_update_interval = weight_update_interval,
           first_factor_update    = first_factor_update,
           first_weight_update    = first_weight_update,
           sample_all_initially   = sample_all_initially,
           n_fss_updates          = n_fss_updates,
           adapt_factors          = adapt_factors,
           initial_widths         = initial_widths,
           initial_factors        = initial_factors)
      }
