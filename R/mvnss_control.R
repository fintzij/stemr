#' Generate a list of settings for multivariate normal slice sampling
#'
#' @param n_updates number of updates per iteration.
#' @param initial_bracket_width initial width of the slice bracket
#' @param bracket_limits limits for the slice bracket width, defaults to
#'   (0,Inf).
#' @param scale_constant constant multiple of the adaptations determined by
#'  \code{scale_cooling}.
#' @param scale_cooling rate at which to cool the adaptation, defaults to 0.5.
#'  Adaptation contributions are governed by a harmonic sequence:
#'  scale_constant/(iteration/step_size+1)^scale_cooling. The
#'  \code{plot_adaptations} function may be used to plot the adaptation factors.
#' @param step_size adaptation increment for each iteration, defaults to 1.
#' @param stop_adaptation iteration at which adaptation should be terminated, 
#'   defaults to 0 for no adaptation. 
#' @param adaptation_offset iteration offset
#' @param nugget nugget for proposal covariance
#' @param nugget_step_size increment for each iteration, defaults to 100/number
#'   of iterations.
#' @param nugget_cooling rate at which to cool the nugget, defaults to 0.9.
#'   Adaptation contributions are governed by a harmonic sequence:
#'   scale_constant/(iteration/step_size+1)^scale_cooling. The
#'   \code{plot_adaptations} function may be used to plot the adaptation
#'   factors.
#'
#' @return list with control settings for multivariate normal slice sampling
#' @export
mvnss_control <-
      function(n_updates = 1,
               bracket_width = 1,
               bracket_limits = c(0,Inf),
               scale_constant = 1,
               scale_cooling = 2/3,
               step_size = 1,
               stop_adaptation = 0,
               adaptation_offset = 0,
               nugget = NULL,
               nugget_cooling = 0.9, 
               nugget_step_size = NULL) {
          
      if(scale_cooling <=0.5 | scale_cooling > 1) {
          warning("The cooling rate must be between 0.5 and 1.")
      }
      
      list(n_updates             = n_updates,
           bracket_width         = bracket_width,
           bracket_limits        = bracket_limits,
           scale_constant        = scale_constant,
           scale_cooling         = scale_cooling,
           step_size             = step_size,
           stop_adaptation       = stop_adaptation,
           adaptation_offset     = adaptation_offset,
           nugget                = nugget,
           nugget_cooling        = nugget_cooling,
           nugget_step_size      = nugget_step_size)
}
