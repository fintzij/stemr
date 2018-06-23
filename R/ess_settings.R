#' Generates a list of settings for sampling the latent LNA paths and
#' time-varying parameters via elliptical slice sampling.
#'
#' @param n_ess_updates number of elliptical slice sampling updates per MCMC
#'   iteration, defaults to 1.
#' @param ess_bracket_width width of the elliptical slice sampling bracket for
#'   LNA path/joint updates, must be in (0,2*pi], and defaults to 2*pi.
#' @param initdist_bracket_width width of the elliptical slice sampling bracket
#'   for initial distribution updates, must be in (0,2*pi], and defaults to
#'   2*pi.
#' @param tparam_bracket_width width of the elliptical slice sampling bracket
#'   for time varying parameter updates, must be in (0,2*pi], and defaults to
#'   2*pi.
#' @param joint_tparam_update either TRUE (default) if time-varying parameter
#'   values should be updated in an ESS step jointly with the LNA path, or FALSE
#'   if time-varying parameters should be updated in their own block.
#' @param joint_initdist_update either TRUE (default) if initial compartment
#'   volumes should be updated in an ESS step jointly with the LNA path, or
#'   FALSE if initial compartment volumes should be updated in their own block.
#'   ess_bracket_step_size+1)^-ess_bracket_cooling.
#' @param ess_warmup ess_warmup ESS updates prior to starting MCMC
#'
#' @return list with settings for elliptical slice sampling
#' @export
ess_settings <-
      function(n_ess_updates = 1,
               ess_bracket_width = 2*pi,
               initdist_bracket_width = 2*pi,
               tparam_bracket_width = 2*pi,
               joint_tparam_update = TRUE,
               joint_initdist_update = TRUE,
               ess_warmup = 50) {
            
            if(ess_bracket_width <= 0 | ess_bracket_width > 2*pi) {
                  stop("The elliptical slice sampling bracket width must be in (0,2*pi].")
            }
            
            if(initdist_bracket_width <= 0 | initdist_bracket_width > 2*pi) {
                  stop("The initial distribution slice sampling bracket width must be in (0,2*pi].")
            }
            
            if(ess_bracket_width <= 0 | ess_bracket_width > 2*pi) {
                  stop("The time varying parameter slice sampling bracket width must be in (0,2*pi].")
            }
            
        return(list(n_ess_updates          = n_ess_updates,
                    ess_bracket_width      = ess_bracket_width,
                    initdist_bracket_width = initdist_bracket_width,
                    tparam_bracket_width   = tparam_bracket_width,
                    tparam_update          = joint_tparam_update,
                    initdist_update        = joint_initdist_update,
                    ess_warmup             = ess_warmup))

}