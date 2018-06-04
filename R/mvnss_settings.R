#' Generate a list of settings for automated factor slice sampling
#'
#' @param n_mvnss_updates number of non-isotropoc hit-and-run updates.
#' @param initial_bracket_width initial width of the slice bracket
#' @param cov_update_interval how often the empirical covariance matrix should be updated
#' @param bracket_limits limits for the slice bracket width, defaults to (0,Inf). 
#'
#' @return list with additional settings for hit and run slice sampling
#' @export
mvnss_settings <-
      function(n_mvnss_updates = 1,
               initial_bracket_width = 1,
               cov_update_interval = 1,
               bracket_limits = c(0,Inf)) {
      
      list(n_mvnss_updates       = n_mvnss_updates,
           initial_bracket_width = initial_bracket_width,
           cov_update_interval   = cov_update_interval,
           bracket_limits        = bracket_limits)
}
