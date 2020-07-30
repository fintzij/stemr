#' Sample a new LNA path via elliptical slice sampling.
#'
#' @inheritParams mvnss_update
#' @param lna_ess_schedule
#'
#' @return list with an updated LNA path along with its stochastic
#'   perturbations, the observed dat log-likelihood, and the lna
#'   log-likelihood, and a record of the number of elliptical slice sampling
#'   proposals
#' @export
lna_update <-
    function(path,
             dat,
             params_cur,
             lna_ess_schedule,
             initdist_objects,
             tparam,
             pathmat_prop,
             censusmat,
             draws_prop,
             ess_draws_prop,
             emitmat,
             flow_matrix,
             stoich_matrix,
             times,
             forcing_inds,
             forcing_tcov_inds,
             forcings_out,
             forcing_transfers,
             param_inds,
             const_inds,
             tcovar_inds,
             initdist_inds,
             param_update_inds,
             census_indices,
             event_inds,
             measproc_indmat,
             svd_d,
             svd_U,
             svd_V,
             proc_pointer,
             set_pars_pointer,
             d_meas_pointer,
             do_prevalence,
             joint_initdist_update,
             return_ess_rec,
             step_size) {
        
        # order of strata updates
        ess_order <- sample.int(length(lna_ess_schedule))
        
        # return ess record?
        if(return_ess_rec) {
            for(s in seq_along(lna_ess_schedule)) {
                reset_vec(lna_ess_schedule[[s]]$steps, 1.0)
                reset_vec(lna_ess_schedule[[s]]$angles, 0.0)
            }
        }
        
        # perform the elliptical slice sampling updates
        for(j in ess_order) {
            for(k in seq_len(lna_ess_schedule[[j]]$n_updates)) {
            
                # sample a new set of stochastic perturbations
                ess_draws_prop[lna_ess_schedule[[j]]$ess_inds,] <- 
                    rnorm(ess_draws_prop[lna_ess_schedule[[j]]$ess_inds,])
                
                # choose a likelihood threshold
                threshold <- path$data_log_lik + log(runif(1))
                
                # initial proposal, which also defines a bracket
                # theta <- runif(1, 0, ess_bracket_width)
                # lower <- theta - ess_bracket_width; upper <- theta
                pos <- runif(1) 
                lower <- -lna_ess_schedule[[j]]$bracket_width[j] * pos
                upper <- lower + lna_ess_schedule[[j]]$bracket_width[j]
                theta <- runif(1, lower, upper)
                    
                # initialize the data log likelihood for the proposed path
                data_log_lik_prop <- NULL
                
                # propose a new initial state
                if(joint_initdist_update) {
                    
                    # indices of strata to sample
                    initdist_codes <- lna_ess_schedule[[j]]$initdist_codes
                    bad_draws      <- vector("logical", length(initdist_codes))
                    
                    # choose an ellipse
                    for(s in initdist_codes) {
                        
                        # if the state is not fixed draw new values
                        if(!initdist_objects[[s]]$fixed) {
                            
                            # draw N(0,1)
                            draw_normals(initdist_objects[[s]]$draws_prop)
                            
                            # compute the linear combination
                            copy_vec(dest = initdist_objects[[s]]$draws_ess, 
                                     orig = 
                                         cos(theta) * initdist_objects[[s]]$draws_cur +
                                         sin(theta) * initdist_objects[[s]]$draws_prop)
                            
                            # map to volumes
                            copy_vec(dest = initdist_objects[[s]]$init_volumes_prop,
                                     orig = c(initdist_objects[[s]]$comp_mean +
                                                  c(initdist_objects[[s]]$comp_sqrt_cov %*%
                                                        initdist_objects[[s]]$draws_ess)))
                            
                            # check boundary conditions
                            bad_draws[s] <- 
                                any(initdist_objects[[s]]$init_volumes_prop < 0 |
                                        initdist_objects[[s]]$init_volumes_prop > 
                                        initdist_objects[[s]]$comp_size)
                        }
                    }
                }
                
                if(joint_initdist_update && any(bad_draws)) {
                    data_log_lik_prop <- -Inf
                    
                } else {
                    
                    # copy the new initial compartment counts
                    if(joint_initdist_update) {
                        insert_initdist(parmat = params_cur,
                                        initdist_objects = initdist_objects[initdist_codes],
                                        prop = TRUE,
                                        rowind = 0,
                                        mcmc_rec = FALSE)
                    }
                    
                    # construct the first proposal
                    copy_2_rows(dest = draws_prop,
                                orig = cos(theta) * path$draws[lna_ess_schedule[[j]]$ess_inds, ] + 
                                    sin(theta) * ess_draws_prop[lna_ess_schedule[[j]]$ess_inds, ],
                                inds = lna_ess_schedule[[j]]$ess_inds - 1)
                    
                    # strata not resampled
                    copy_2_rows(dest = draws_prop,
                                orig = path$draws[lna_ess_schedule[[j]]$complementary_inds,],
                                inds = lna_ess_schedule[[j]]$complementary_inds - 1)
                    
                    try({
                    
                        # map draws onto a latent path
                        map_draws_2_lna(
                            pathmat           = pathmat_prop,
                            draws             = draws_prop,
                            lna_times         = times,
                            lna_pars          = params_cur,
                            lna_param_vec     = param_vec,
                            lna_param_inds    = param_inds,
                            lna_tcovar_inds   = tcovar_inds,
                            init_start        = initdist_inds[1],
                            param_update_inds = param_update_inds,
                            stoich_matrix     = stoich_matrix,
                            forcing_inds      = forcing_inds,
                            forcing_tcov_inds = forcing_tcov_inds,
                            forcings_out      = forcings_out,
                            forcing_transfers = forcing_transfers,
                            svd_d             = svd_d,
                            svd_U             = svd_U,
                            svd_V             = svd_V,
                            lna_pointer       = proc_pointer,
                            set_pars_pointer  = set_pars_pointer,
                            step_size         = step_size
                        ) 
                        
                        census_latent_path(
                            path                = pathmat_prop,
                            census_path         = censusmat,
                            census_inds         = census_indices,
                            event_inds          = event_inds,
                            flow_matrix         = flow_matrix,
                            do_prevalence       = do_prevalence,
                            parmat              = params_cur,
                            initdist_inds       = initdist_inds,
                            forcing_inds        = forcing_inds,
                            forcing_tcov_inds   = forcing_tcov_inds,
                            forcings_out        = forcings_out,
                            forcing_transfers   = forcing_transfers
                        )
                        
                        # evaluate the density of the incidence counts
                        evaluate_d_measure_LNA(
                            emitmat           = emitmat,
                            obsmat            = dat,
                            censusmat         = censusmat,
                            measproc_indmat   = measproc_indmat,
                            parameters        = params_cur,
                            param_inds        = param_inds,
                            const_inds        = const_inds,
                            tcovar_inds       = tcovar_inds,
                            param_update_inds = param_update_inds,
                            census_indices    = census_indices,
                            param_vec         = param_vec,
                            d_meas_ptr        = d_meas_pointer)
                        
                        # compute the data log likelihood
                        data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                        if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
                    }, silent = TRUE)
                    
                    # if proposal failed data_log_lik_prop is -Inf
                    if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
                }
                
                # continue proposing if not accepted
                while((upper - lower) > sqrt(.Machine$double.eps) && (data_log_lik_prop < threshold)) {
                    
                    # increment the number of ESS steps
                    if(return_ess_rec) 
                        increment_elem(lna_ess_schedule[[j]]$steps, k)
                    
                    # shrink the bracket
                    if(theta < 0) {
                        lower <- theta
                    } else {
                        upper <- theta
                    }
                    
                    # sample a new point
                    theta <- runif(1, lower, upper)
                    
                    # construct the next initial distribution proposal
                    if(joint_initdist_update) {
                        for(s in initdist_codes) {
                            
                            # if the state is not fixed draw new values
                            if(!initdist_objects[[s]]$fixed) {
                                
                                # compute the linear combination
                                copy_vec(dest = initdist_objects[[s]]$draws_ess, 
                                         orig = 
                                             cos(theta) * initdist_objects[[s]]$draws_cur +
                                             sin(theta) * initdist_objects[[s]]$draws_prop)
                                
                                # map to volumes
                                copy_vec(dest = initdist_objects[[s]]$init_volumes_prop,
                                         orig = c(initdist_objects[[s]]$comp_mean +
                                                      c(initdist_objects[[s]]$comp_sqrt_cov %*%
                                                            initdist_objects[[s]]$draws_ess)))
                                
                                # check boundary conditions
                                bad_draws[s] <- 
                                    any(initdist_objects[[s]]$init_volumes_prop < 0 |
                                            initdist_objects[[s]]$init_volumes_prop > 
                                            initdist_objects[[s]]$comp_size)
                            }
                        }
                    }
                    
                    if(joint_initdist_update && any(bad_draws)) {
                        data_log_lik_prop <- -Inf
                        
                    } else {
                        
                        # copy the new initial compartment counts
                        if(joint_initdist_update) {
                            insert_initdist(parmat = params_cur,
                                            initdist_objects = initdist_objects[initdist_codes],
                                            prop = TRUE,
                                            rowind = 0,
                                            mcmc_rec = FALSE)
                        }
                        
                        # construct the first proposal
                        copy_2_rows(dest = draws_prop,
                                    orig = cos(theta) * path$draws[lna_ess_schedule[[j]]$ess_inds, ] + 
                                        sin(theta) * ess_draws_prop[lna_ess_schedule[[j]]$ess_inds, ],
                                    inds = lna_ess_schedule[[j]]$ess_inds - 1)
                        
                        # strata not resamples
                        copy_2_rows(dest = draws_prop,
                                    orig = path$draws[lna_ess_schedule[[j]]$complementary_inds,],
                                    inds = lna_ess_schedule[[j]]$complementary_inds - 1)
                        
                        try({
                            # map draws onto a latent path
                            map_draws_2_lna(
                                pathmat           = pathmat_prop,
                                draws             = draws_prop,
                                lna_times         = census_times,
                                lna_pars          = params_cur,
                                lna_param_vec     = param_vec,
                                lna_param_inds    = param_inds,
                                lna_tcovar_inds   = tcovar_inds,
                                init_start        = initdist_inds[1],
                                param_update_inds = param_update_inds,
                                stoich_matrix     = stoich_matrix,
                                forcing_inds      = forcing_inds,
                                forcing_tcov_inds = forcing_tcov_inds,
                                forcings_out      = forcings_out,
                                forcing_transfers = forcing_transfers,
                                svd_d             = svd_d,
                                svd_U             = svd_U,
                                svd_V             = svd_V,
                                lna_pointer       = proc_pointer,
                                set_pars_pointer  = set_pars_pointer,
                                step_size         = step_size
                            ) 
                            
                            census_latent_path(
                                path                = pathmat_prop,
                                census_path         = censusmat,
                                census_inds         = census_indices,
                                event_inds          = event_inds,
                                flow_matrix         = flow_matrix,
                                do_prevalence       = do_prevalence,
                                parmat              = params_cur,
                                initdist_inds       = initdist_inds,
                                forcing_inds        = forcing_inds,
                                forcing_tcov_inds   = forcing_tcov_inds,
                                forcings_out        = forcings_out,
                                forcing_transfers   = forcing_transfers
                            )
                            
                            # evaluate the density of the incidence counts
                            evaluate_d_measure_LNA(
                                emitmat           = emitmat,
                                obsmat            = dat,
                                censusmat         = censusmat,
                                measproc_indmat   = measproc_indmat,
                                parameters        = params_cur,
                                param_inds        = param_inds,
                                const_inds        = const_inds,
                                tcovar_inds       = tcovar_inds,
                                param_update_inds = param_update_inds,
                                census_indices    = census_indices,
                                param_vec         = param_vec,
                                d_meas_ptr        = d_meas_pointer)
                            
                            # compute the data log likelihood
                            data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                            if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
                        }, silent = TRUE)
                        
                        # if proposal failed data_log_lik_prop is -Inf
                        if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
                    }
                }
                
                # if the bracket width is not equal to zero, update the draws, path, and dat log likelihood
                if((upper - lower) > sqrt(.Machine$double.eps)) {
                    
                    # transfer the new initial volumes and draws (volumes already in parameter matrix)
                    if(joint_initdist_update) {
                        for(s in initdist_codes) {
                            if(!initdist_objects[[s]]$fixed) {
                                
                                # copy the N(0,1) draws
                                copy_vec(dest = initdist_objects[[s]]$draws_cur,
                                         orig = initdist_objects[[s]]$draws_ess)
                                
                                # copy the initial compartment volumes
                                copy_vec(dest = initdist_objects[[s]]$init_volumes,
                                         orig = initdist_objects[[s]]$init_volumes_prop)
                            }
                        }
                    }
                    
                    # transfer the new path and residual path into the path list
                    copy_2_rows(dest = path$draws,
                                orig = draws_prop[lna_ess_schedule[[j]]$ess_inds,],
                                inds = lna_ess_schedule[[j]]$ess_inds-1)
                    
                    # copy the LNA path and the dat log likelihood
                    copy_mat(dest = path$latent_path, orig = pathmat_prop)
                    copy_vec(dest = path$data_log_lik, orig = data_log_lik_prop)
                    
                    # record the final angle
                    if(return_ess_rec) {
                        insert_elem(dest = lna_ess_schedule[[j]]$angles,
                                    elem = theta,
                                    ind  = k-1)
                    }
                    
                } else {
                    
                    # insert the original compartment counts back into the parameter matrix
                    if(joint_initdist_update) {
                        insert_initdist(parmat = params_cur,
                                        initdist_objects = initdist_objects[initdist_codes],
                                        prop = FALSE,
                                        rowind = 0,
                                        mcmc_rec = FALSE)
                    }
                }   
            }
        }
    }
