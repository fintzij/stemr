#' Generates a list of settings for sampling the latent LNA paths via elliptical
#' slice sampling.
#'
#' @param n_ess_updates number of elliptical slice sampling updates per MCMC
#'   iteration, defaults to 1.
#' @param tparam_update either "joint" (default) if time-varying parameter
#'   values should be updated jointly with the LNA path, or "block" if
#'   time-varying parameters should be updated in their own block.
#' @param ess_warmup ess_warmup ESS updates prior to starting MCMC
#'
#' @return list with settings for elliptical slice sampling
#' @export
ess_settings <-
      function(n_ess_updates = 1,
               tparam_update = "joint",
               ess_warmup = 50) {
            
        return(list(n_ess_updates = n_ess_updates,
                    tparam_update = tparam_update,
                    ess_warmup    = ess_warmup))

}