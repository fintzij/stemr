#' Save an MCMC sample
#'
#' @param mcmc_samples list with objects for recording MCMC samples
#' @param rec_ind C++ record index
#' @param path latent path
#' @param param_blocks list of parameter blocks
#' @param initdist_objects list of initial distribution objects
#' @param tparam list with time-varying parameters
#' @param method method for approximating the latent epidemic, either "lna" or
#'   "ode" 
#'
#' @return inserts all parameters into their respective MCMC record
#' @export
save_mcmc_sample = 
    function(mcmc_samples,
             rec_ind,
             path,
             parmat,
             param_blocks, 
             initdist_objects,
             tparam, 
             tparam_inds,
             method) {
        
        # record parameters
        insert_params(parmat       = mcmc_samples$parameter_samples_nat,
                      param_blocks = param_blocks,
                      nat          = TRUE,
                      prop         = FALSE,
                      rowind       = rec_ind)
        
        insert_params(parmat       = mcmc_samples$parameter_samples_est,
                      param_blocks = param_blocks,
                      nat          = FALSE,
                      prop         = FALSE,
                      rowind       = rec_ind)
        
        insert_elem(dest = mcmc_samples$params_log_prior,
                    elem = sum(sapply(param_blocks, "[[", "log_pd")),
                    ind  = rec_ind)
        
        # copy latent path
        mat_2_arr(dest = mcmc_samples$latent_paths,
                  orig = path$latent_path,
                  ind  = rec_ind)
        
        insert_elem(dest = mcmc_samples$data_log_lik,
                    elem = path$data_log_lik,
                    ind  = rec_ind)
        
        if(method == "lna") {
            mat_2_arr(dest = mcmc_samples$lna_draws,
                      orig = path$draws,
                      ind = rec_ind)
            
            insert_elem(dest = mcmc_samples$lna_log_lik,
                        elem = sum(dnorm(path$draws, log = TRUE)),
                        ind  = rec_ind)
        }
        
        # record initial compartment counts
        if(!is.null(mcmc_samples$initdist_samples)) {
            insert_initdist(mcmc_samples$initdist_samples,
                            initdist_objects,
                            prop = FALSE, 
                            rowind = rec_ind,
                            mcmc_rec = TRUE)
            
            insert_elem(dest = mcmc_samples$initdist_log_lik,
                        elem = sum(dnorm(sapply(initdist_objects, "[[", "draws_cur"), log = T)),
                        ind  = rec_ind)
        }
        
        # record time-varying parameters
        if(!is.null(tparam_inds)) {
            mat_2_arr(dest = mcmc_samples$tparam_samples,
                      orig = parmat[, tparam_inds + 1, drop=FALSE],
                      ind  = rec_ind)
            
            insert_elem(dest = mcmc_samples$tparam_log_lik,
                        elem = sum(dnorm(sapply(tparam, "[[", "draws_cur"), log = T)),
                        ind  = rec_ind)
        }
        
        # increment rec_ind
        increment_vec(rec_ind, 1)
    }
