#' Generates a list of settings for sampling the latent LNA paths and
#' time-varying parameters via elliptical slice sampling.
#'
#' @param n_updates number of elliptical slice sampling updates per MCMC
#'   iteration, defaults to 1.
#' @param bracket_width width of the elliptical slice sampling bracket, must
#'   be in (0,2*pi], and default to 2*pi. 
#' @param bracket_update_iter iteration at which the widths of the elliptical
#'   slice sampling brackets should be shrunk, defaults to Inf and the bracket
#'   widths are kept constant.
#' @param bracket_scaling Scaling factors for elliptical slice sampling bracket
#'   widths. If brackets are to be shrunk, the new width is set to the minimum
#'   of 2*pi or the scaling multiplied by the standard deviation of the previous
#'   ESS angles. The scalings default to the full width at one tenth maximum for
#'   a gaussian.
#' @param joint_strata_update should all strata be updated jointly? Defaults to
#'   FALSE.
#' @param joint_initdist_update should the initial states be updated jointly
#'   with the lna path? Defaults to TRUE, in which case initial conditions for
#'   each stratum are still paired with the LNA path for that stratum.
#' @param ess_warmup number of ESS warmup iterations if approximate LNA
#'   initialization
#'
#' @return list with settings for elliptical slice sampling
#' @export
lna_control <-
      function(n_updates = 1,
               bracket_width = 2 * pi,
               bracket_update_iter = Inf,
               bracket_scaling = 2 * sqrt(2 * log(10)),
               joint_strata_update = FALSE,
               joint_initdist_update = TRUE,
               ess_warmup = 100) {
          
            if (any(bracket_width <= 0 | bracket_width > 2 * pi)) {
                  stop("The elliptical slice sampling bracket width must be in (0,2*pi].")
            }
            
            return(
                  list(n_updates             = n_updates,
                       bracket_width         = bracket_width,
                       bracket_update_iter   = bracket_update_iter,
                       bracket_scaling       = bracket_scaling,
                       joint_strata_update   = joint_strata_update,
                       joint_initdist_update = joint_initdist_update,
                       ess_warmup            = ess_warmup
                  )
            )
      }