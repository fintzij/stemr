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
#'   sqrt(lambda_i)/sum(sqrt(lambda_i))$, where $lambda_i$ are the eigen values
#'   of the empirical covariance matrix.
#' @param first_factor_update iteration at which the first factor update should
#'   occur
#' @param first_prob_update iteration at which the first update to the slice
#'   direction probabilities should occur.
#' @param initial_slice_probs vector of initial slice probabilities
#' @param n_afss_updates number of factor slice sampling directions to update
#' @param target_prop_totsd Between 0 and 1. If supplied, the number of factor
#'   directions sampled after the first probability update set to be the maximum
#'   of the requested number of directions and number of factors required to
#'   explain the specified proportion of the total standard deviation in the
#'   posterior. If not equal to the number of parameters, a hit and run update
#'   is carried out to ensure ergodicity.
#' @param sample_all_initially should all factor directions be sampled until the
#'   first slice probability update? defaults to TRUE
#' @param harss_prob probability of a hit and run update at each iteration
#' @param harss_warmup should a hit and run update be performed along with the
#'   elliptical slice sampling warm up? defaults to TRUE
#' @param initial_widths vector of initial slice widths, defaults to a vector of
#' @param afss_slice_ratio target ratio of expansions / (expansions +
#'   contractions) ones.
#' @return list with additional settings for automated factor slice sampling
#' @export
afss_settings <-
      function(factor_update_interval = NULL,
               prob_update_interval = NULL,
               first_factor_update = 100,
               first_prob_update = 100,
               initial_slice_probs = NULL,
               initial_widths = NULL,
               n_afss_updates = NULL,
               sample_all_initially = TRUE,
               afss_slice_ratio = 0.5,
               target_prop_totsd = NULL,
               harss_prob = 0.05,
               harss_warmup = TRUE) {
               
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
            
      if(!is.null(target_prop_totsd) && !(target_prop_totsd >= 0 & target_prop_totsd < 1)) {
            stop("The proportion of the total standard deviation explained by the afss proposals must be between 0 and 1.")
      }
      
      if(is.null(prob_update_interval)) prob_update_interval <- factor_update_interval
      
      list(factor_update_interval  = factor_update_interval,
           prob_update_interval    = prob_update_interval,
           first_factor_update     = first_factor_update,
           first_prob_update       = first_prob_update,
           initial_slice_probs     = initial_slice_probs,
           initial_widths          = initial_widths,
           n_afss_updates          = n_afss_updates,
           sample_all_initially    = sample_all_initially,
           afss_slice_ratio        = afss_slice_ratio,
           target_prop_totsd       = target_prop_totsd,
           harss_prob              = harss_prob,
           harss_warmup            = harss_warmup)
      }
