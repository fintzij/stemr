#' Sample a new LNA path via elliptical slice sampling.
#'
#' @inheritParams mvnss_update
#'
#' @return updates the initial conditions
#' @export
initdist_update <-
    function(path,
             dat,
             iter,
             parmat,
             initdist_objects,
             initdist_ess_control,
             tparam,
             pathmat_prop,
             censusmat,
             draws_prop,
             ess_draws_prop,
             emitmat,
             flow_matrix,
             stoich_matrix,
             census_times,
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
             svd_d = NULL,
             svd_U = NULL,
             svd_V = NULL,
             proc_pointer,
             set_pars_pointer,
             d_meas_pointer,
             do_prevalence,
             step_size) {
        
    # reset ESS steps and angles
    reset_vec(initdist_ess_control$steps, 1.0)
    reset_vec(initdist_ess_control$angles, 0.0)
    
    # perform the elliptical slice sampling updates
    for(k in seq_len(initdist_ess_control$n_updates)) {
        
        # choose a likelihood threshold
        threshold <- path$data_log_lik + log(runif(1))
        
        # initialize the data log likelihood for the proposed path
        data_log_lik_prop <- NULL
        
        # initial proposal, which also defines a bracket
        pos <- runif(1) 
        lower <- -initdist_ess_control$bracket_width * pos
        upper <- lower + initdist_ess_control$bracket_width
        theta <- runif(1, lower, upper)
        
        # vector of logicals for whether boundary conditions are respected
        bad_draws <- vector("logical", length(initdist_objects))
        
        # choose an ellipse
        for(s in seq_along(initdist_objects)) {
            
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
            
        if(any(bad_draws)) {
            data_log_lik_prop <- -Inf
            
        } else {
            
            # copy the new initial compartment counts
            insert_initdist(parmat = parmat,
                            initdist_objects = initdist_objects,
                            prop = TRUE,
                            rowind = 0,
                            mcmc_rec = FALSE)
            
            if(!is.null(tparam)) {
                for(p in seq_along(tparam)) {
                    if(tparam[[p]]$init_dep) {
                        
                        insert_tparam(
                            tcovar = parmat,
                            values = 
                                tparam[[p]]$draws2par(
                                    parameters = parmat[1,],
                                    draws = tparam[[p]]$draws_cur),
                            col_ind = tparam[[p]]$col_ind,
                            tpar_inds = tparam[[p]]$tpar_inds_Cpp)    
                    }
                }
            }
            
            # map the perturbations to a latent path
            try({
                if(is.null(svd_d)) {
                    
                    map_pars_2_ode(
                        pathmat           = pathmat_prop,
                        ode_times         = census_times,
                        ode_pars          = parmat,
                        ode_param_vec     = param_vec,
                        ode_param_inds    = param_inds,
                        ode_tcovar_inds   = tcovar_inds,
                        init_start        = initdist_inds[1],
                        param_update_inds = param_update_inds,
                        stoich_matrix     = stoich_matrix,
                        forcing_inds      = forcing_inds,
                        forcing_tcov_inds = forcing_tcov_inds,
                        forcings_out      = forcings_out,
                        forcing_transfers = forcing_transfers,
                        ode_pointer       = proc_pointer,
                        set_pars_pointer  = set_pars_pointer,
                        step_size         = step_size
                    ) 
                    
                } else {
                    map_draws_2_lna(
                        pathmat           = pathmat_prop,
                        draws             = path$draws,
                        lna_times         = census_times,
                        lna_pars          = parmat,
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
                }
                
                census_latent_path(
                    path                = pathmat_prop,
                    census_path         = censusmat,
                    census_inds         = census_indices,
                    event_inds          = event_inds,
                    flow_matrix         = flow_matrix,
                    do_prevalence       = do_prevalence,
                    parmat              = parmat,
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
                    parameters        = parmat,
                    param_inds        = param_inds,
                    const_inds        = const_inds,
                    tcovar_inds       = tcovar_inds,
                    param_update_inds = param_update_inds,
                    census_indices    = census_indices,
                    param_vec         = param_vec,
                    d_meas_ptr        = d_meas_pointer)
                
                # compute the data log likelihood
                data_log_lik_prop <- sum(emitmat[, -1][measproc_indmat])
                
                if (is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
            }, silent = TRUE)
            
            # if proposal failed data_log_lik_prop is -Inf
            if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
        }
            
        # continue proposing if not accepted
        while((upper - lower) > sqrt(.Machine$double.eps) && 
              (data_log_lik_prop < threshold)) {
            
            # increment the number of ESS steps
            increment_elem(initdist_ess_control$steps, k-1)
            
            # shrink the bracket
            if(theta < 0) {
                lower <- theta
            } else {
                upper <- theta
            }
            
            # sample a new point
            theta <- runif(1, lower, upper)
            
            # construct the next initial distribution proposal
            for(s in seq_along(initdist_objects)) {
                
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
            
            if(any(bad_draws)) {
                data_log_lik_prop <- -Inf
                
            } else {
                
                # copy the new initial compartment counts
                insert_initdist(parmat = parmat,
                                initdist_objects = initdist_objects,
                                prop = TRUE,
                                rowind = 0,
                                mcmc_rec = FALSE)
                
                if(!is.null(tparam)) {
                    for(p in seq_along(tparam)) {
                        if(tparam[[p]]$init_dep) {
                            
                            insert_tparam(
                                tcovar = parmat,
                                values = 
                                    tparam[[p]]$draws2par(
                                        parameters = parmat[1,],
                                        draws = tparam[[p]]$draws_cur),
                                col_ind = tparam[[p]]$col_ind,
                                tpar_inds = tparam[[p]]$tpar_inds_Cpp)    
                        }
                    }
                }
                
                # map the perturbations to a latent path
                try({
                    if(is.null(svd_d)) {
                        
                        map_pars_2_ode(
                            pathmat           = pathmat_prop,
                            ode_times         = census_times,
                            ode_pars          = parmat,
                            ode_param_vec     = param_vec,
                            ode_param_inds    = param_inds,
                            ode_tcovar_inds   = tcovar_inds,
                            init_start        = initdist_inds[1],
                            param_update_inds = param_update_inds,
                            stoich_matrix     = stoich_matrix,
                            forcing_inds      = forcing_inds,
                            forcing_tcov_inds = forcing_tcov_inds,
                            forcings_out      = forcings_out,
                            forcing_transfers = forcing_transfers,
                            ode_pointer       = proc_pointer,
                            set_pars_pointer  = set_pars_pointer,
                            step_size         = step_size
                        ) 
                        
                    } else {
                        map_draws_2_lna(
                            pathmat           = pathmat_prop,
                            draws             = path$draws,
                            lna_times         = census_times,
                            lna_pars          = parmat,
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
                    }
                    
                    census_latent_path(
                        path                = pathmat_prop,
                        census_path         = censusmat,
                        census_inds         = census_indices,
                        event_inds          = event_inds,
                        flow_matrix         = flow_matrix,
                        do_prevalence       = do_prevalence,
                        parmat              = parmat,
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
                        parameters        = parmat,
                        param_inds        = param_inds,
                        const_inds        = const_inds,
                        tcovar_inds       = tcovar_inds,
                        param_update_inds = param_update_inds,
                        census_indices    = census_indices,
                        param_vec         = param_vec,
                        d_meas_ptr        = d_meas_pointer)
                    
                    # compute the data log likelihood
                    data_log_lik_prop <- sum(emitmat[, -1][measproc_indmat])
                    
                    if (is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
                }, silent = TRUE)
                
                # if proposal failed data_log_lik_prop is -Inf
                if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
            }
        }
                
        # if the bracket width is not equal to zero, update the draws, path, and dat log likelihood
        if((upper - lower) > sqrt(.Machine$double.eps)) {
            
            # transfer the new initial volumes and draws (volumes already in parameter matrix)
            for(s in seq_along(initdist_objects)) {
                if(!initdist_objects[[s]]$fixed) {
                    
                    # copy the N(0,1) draws
                    copy_vec(dest = initdist_objects[[s]]$draws_cur,
                             orig = initdist_objects[[s]]$draws_ess)
                    
                    # copy the initial compartment volumes
                    copy_vec(dest = initdist_objects[[s]]$init_volumes,
                             orig = initdist_objects[[s]]$init_volumes_prop)
                }
            }
            
            # copy the LNA path and the dat log likelihood
            copy_vec(dest = path$data_log_lik, orig = data_log_lik_prop)
            
            if(k == initdist_ess_control$n_updates) {
                copy_mat(dest = path$latent_path, orig = pathmat_prop)
            }
            
            # copy time-varying parameters
            if(!is.null(tparam)) {
                for(p in seq_along(tparam)) {
                    if(tparam[[p]]$init_dep) {
                        copy_vec(dest = tparam[[p]]$tpar_cur,
                                 orig = parmat[,tparam[[p]]$col_ind + 1])    
                    }
                }
            }
            
            # record the final angle
            insert_elem(dest = initdist_ess_control$angles,
                        elem = theta,
                        ind  = k-1)
            
        } else {
            
            # insert the original compartment counts back into the parameter matrix
            insert_initdist(parmat = parmat,
                            initdist_objects = initdist_objects,
                            prop = FALSE,
                            rowind = 0,
                            mcmc_rec = FALSE)
            
            # recover the original time-varying parameter values
            if(!is.null(tparam)) {
                for(p in seq_along(tparam)) {
                    if(tparam[[p]]$init_dep) {
                        vec_2_mat(dest = parmat,
                                  orig = tparam[[p]]$tpar_cur,
                                  ind = tparam[[p]]$col_ind)
                    }
                }
            }
        }
    }
    
    # update the ESS bracket
    if(iter != 0 && iter <= initdist_ess_control$bracket_update_iter) {
        
        # angle residual
        copy_vec(
            dest = initdist_ess_control$angle_resid,
            orig = mean(initdist_ess_control$angles) - 
                initdist_ess_control$angle_mean) 
        
        # angle variance
        copy_vec(
            dest = initdist_ess_control$angle_var,
            orig = initdist_ess_control$angle_resid^2 / iter + 
                initdist_ess_control$angle_var * (iter-1) / iter)
        
        # angle mean
        copy_vec(
            dest = initdist_ess_control$angle_mean,
            orig = mean(initdist_ess_control$angles) / iter + 
                initdist_ess_control$angle_mean * (iter-1) / iter)
        
        # set the new angle bracket
        if(iter == initdist_ess_control$bracket_update_iter) {
            copy_vec(
                dest = initdist_ess_control$bracket_width,
                orig = pmin(initdist_ess_control$bracket_scaling * 
                                sqrt(initdist_ess_control$angle_var),
                            2*pi))    
        }
    }
}
    
