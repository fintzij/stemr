#' Generates a list of settings for sampling the latent LNA paths and
#' time-varying parameters via elliptical slice sampling.
#'
#' @param n_updates number of elliptical slice sampling updates per MCMC
#'   iteration, defaults to 1.
#' @param bracket_widths width of the elliptical slice sampling brackets, must
#'   be in (0,2*pi], and default to 2*pi.
#' @param bracket_update_iter iteration at which the widths of the elliptical
#'   slice sampling brackets should be shrunk, defaults to 0 and the bracket
#'   widths are kept constant.
#' @param bracket_scaling Scaling factors for elliptical slice sampling bracket
#'   widths. If brackets are to be shrunk, the new width is set to the minimum
#'   of 2*pi or the scaling multiplied by the standard deviation of the previous
#'   ESS angles. The scalings default to the full width at one tenth maximum for
#'   a gaussian.
#' @param joint_tparam_update should all time-varying parameters be updated
#'   jointly? Defaults to FALSE.
#' @param ess_warmup ess_warmup ESS updates prior to starting MCMC
#'
#' @return list with settings for elliptical slice sampling
#' @export
tpar_control <-
      function(n_updates = 1,
               bracket_widths = 2 * pi,
               bracket_update_iter = 0,
               bracket_scaling = 2 * sqrt(2 * log(10)),
               joint_tparam_update = FALSE,
               ess_warmup = 50) {
          
            if (any(bracket_widths <= 0 | bracket_widths > 2 * pi)) {
                  stop("The elliptical slice sampling bracket width must be in (0,2*pi].")
            }
            
            return(
                  list(n_updates           = n_updates,
                       bracket_widths      = bracket_widths,
                       bracket_update_iter = bracket_update_iter,
                       bracket_scaling     = bracket_scaling,
                       joint_tparam_update = joint_tparam_update,
                       ess_warmup          = ess_warmup
                  )
            )
      }