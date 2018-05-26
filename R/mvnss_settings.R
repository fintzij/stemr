#' Generate a list of settings for automated factor slice sampling
#'
#' @param n_mvnss_updates number of non-isotropoc hit-and-run updates.
#' @param cov_update_interval how often the empirical covariance matrix should be updated
#'
#' @return list with additional settings for hit and run slice sampling
#' @export
mvnss_settings <-
      function(n_mvnss_updates = 1,
               initial_bracket_width = 1,
               cov_update_interval = 1) {
      
      list(n_mvnss_updates         = n_mvnss_updates,
           initial_bracket_width   = initial_bracket_width,
           cov_update_interval     = cov_update_interval)
}
