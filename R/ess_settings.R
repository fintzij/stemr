#' Generates a list of settings for sampling the latent LNA paths and
#' time-varying parameters via elliptical slice sampling.
#'
#' @param n_ess_updates number of elliptical slice sampling updates per MCMC
#'   iteration, defaults to 1.
#' @param joint_tparam_update either TRUE (default) if time-varying parameter values
#'   should be updated in an ESS step jointly with the LNA path, or FALSE if
#'   time-varying parameters should be updated in their own block.
#' @param joint_initdist_update either TRUE (default) if initial compartment volumes
#'   should be updated in an ESS step jointly with the LNA path, or FALSE if
#'   initial compartment volumes should be updated in their own block.
#' @param ess_warmup ess_warmup ESS updates prior to starting MCMC
#'
#' @return list with settings for elliptical slice sampling
#' @export
ess_settings <-
      function(n_ess_updates = 1,
               joint_tparam_update = TRUE,
               joint_initdist_update = TRUE,
               ess_warmup = 50) {
            
        return(list(n_ess_updates = n_ess_updates,
                    tparam_update = joint_tparam_update,
                    initdist_update = joint_initdist_update,
                    ess_warmup    = ess_warmup))

}