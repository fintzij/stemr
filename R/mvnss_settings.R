#' Generate a list of settings for automated factor slice sampling
#'
#' @param n_mvnss_updates number of non-isotropoc hit-and-run updates.
#' @param initial_bracket_width initial width of the slice bracket
#' @param cov_update_interval how often the empirical covariance matrix should
#'   be updated
#' @param bracket_limits limits for the slice bracket width, defaults to
#'   (0,Inf).
#' @param nugget_cooling rate at which to cool the nugget, defaults to 0.9.
#'   Adaptation contributions are governed by a harmonic sequence:
#'   scale_constant/(iteration/step_size+1)^scale_cooling. The
#'   \code{plot_adaptations} function may be used to plot the adaptation
#'   factors.
#' @param nugget_step_size increment for each iteration, defaults to 100/number
#'   of iterations.
#'
#' @return list with additional settings for hit and run slice sampling
#' @export
mvnss_settings <-
      function(n_mvnss_updates = 1,
               initial_bracket_width = 1,
               cov_update_interval = 1,
               bracket_limits = c(0,Inf),
               nugget_cooling = 0.9, 
               nugget_step_size = NULL) {
      
      list(n_mvnss_updates       = n_mvnss_updates,
           initial_bracket_width = initial_bracket_width,
           cov_update_interval   = cov_update_interval,
           bracket_limits        = bracket_limits,
           nugget_cooling        = nugget_cooling,
           nugget_step_size      = nugget_step_size)
}
