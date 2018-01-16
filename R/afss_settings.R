#' Generate a list of settings for automated factor slice sampling
#'
#' @param factor_update_interval number of iterations between factor and
#'   interval updates (defaults to an update every iteration). May also be
#'   supplied as a function that takes the current iteration as an input and
#'   returns the next iteration at which the factors should be updated. Must
#'   have a default set so that the first update can be determined. e.g.,
#'   \code{function(x = 10) 2*x} will first update the factors at iteration 10
#'   and will subsequently update the factors at iterations 20, 40, 80, ...
#' @param prob_update_interval number of iterations between updates to the
#'   factor sampling probabilities, defaults to the factor update interval. May
#'   also be supplied as a function that takes the current iteration as an input
#'   and returns the next iteration at which the factors should be updated. Must
#'   have a default set so that the first update can be determined. The factor
#'   sampling probabilities are $min(nugget,
#'   sqrt(lambda_i)/sum(sqrt(lambda_i))$, where $lambda_i$ are the singular
#'   values of the empirical covariance matrix.
#' @param initial_widths vector of initial slice widths
#' @param initial_factors matrix of initial factor loadings
#' @param first_factor_update iteration at which the first factor update should
#'   occur
#' @param first_prob_update iteration at which the first update to the slice
#'   direction probabilities should occur.
#' @param min_afss_updates minimum number of afss updates per iteration,
#'   defaults to 1.
#' @param initial_slice_probs initial slice direction probabilities
#' @param use_cor should the slice directions be computed as singular vectors of
#'   the correlation matrix (as opposed to the covariance)? defaults to TRUE.
#'
#' @return list with additional settings for automated factor slice sampling
#' @export
afss_settings <-
      function(factor_update_interval = NULL,
               prob_update_interval = NULL,
               first_factor_update = 100,
               first_prob_update = 100,
               min_afss_updates = NULL,
               initial_widths = NULL,
               initial_factors = NULL,
               initial_slice_probs = NULL,
               use_cor = TRUE) {
         
      # test if the factor and prob update interval functions have defaults if they are specified via functions   
      if(is.function(factor_update_interval)) {
            fui <- NULL
            try({fui <- factor_update_interval()}, silent = T)
            if(is.null(fui)) stop("If the factor update interval is supplied as a function, it must have a default value for the first interval.")
      }
      
      if(is.function(prob_update_interval)) {
            fui <- NULL
            try({fui <- prob_update_interval()}, silent = T)
            if(is.null(fui)) stop("If the probability update interval is supplied as a function, it must have a default value for the first interval.")
      }
      
      if(is.null(prob_update_interval)) prob_update_interval <- factor_update_interval
      
      list(factor_update_interval = factor_update_interval,
           prob_update_interval   = prob_update_interval,
           first_factor_update    = first_factor_update,
           first_prob_update      = first_prob_update,
           min_afss_updates       = min_afss_updates,
           initial_widths         = initial_widths,
           initial_factors        = initial_factors,
           initial_slice_probs    = initial_slice_probs,
           use_cor                = use_cor)
      }
