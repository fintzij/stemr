#' Sample a new LNA path via elliptical slice sampling.
#'
#' @inheritParams mvnss_update
#'
#' @return updates the initial conditions
#' @export
tparam_update <-
    function(path,
             dat,
             iter,
             parmat,
             tparam,
             tparam_ess_control,
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
             svd_d = NULL,
             svd_U = NULL,
             svd_V = NULL,
             proc_pointer,
             set_pars_pointer,
             d_meas_pointer,
             do_prevalence,
             step_size) {
        
    # reset ESS steps and angles
    for(p in seq_along(tparam)) {
        reset_vec(tparam[[p]]$steps, 1.0)
        reset_vec(tparam[[p]]$angles, 0.0)
    }
        
    # perform the elliptical slice sampling updates
    for(k in seq_len(tparam_ess_control$n_updates)) {
        
        # order of tparam updates
        ess_order <- sample.int(length(tparam))
        
        for(p in ess_order) {
    
            # choose a likelihood threshold
            threshold <- path$data_log_lik + log(runif(1))
            
            # initialize the data log likelihood for the proposed path
            data_log_lik_prop <- NULL
            
            # initial proposal, which also defines a bracket
            pos <- runif(1) 
            lower <- -tparam[[p]]$bracket_width * pos
            upper <- lower + tparam[[p]]$bracket_width
            theta <- runif(1, lower, upper)
        
            # construct the first proposal
            draw_normals(tparam[[p]]$draws_prop)
            
            # compute the proposal
            copy_vec(
                dest = tparam[[p]]$draws_ess,
                orig = cos(theta) * tparam[[p]]$draws_cur + 
                       sin(theta) * tparam[[p]]$draws_prop
            )
            
            # insert time-varying parameters
            insert_tparam(
                tcovar = parmat,
                values = 
                    tparam[[p]]$draws2par(
                        parameters = parmat[1,],
                        draws = tparam[[p]]$draws_ess),
                col_ind = tparam[[p]]$col_ind,
                tpar_inds = tparam[[p]]$tpar_inds_Cpp)    
            
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
        
            # continue proposing if not accepted
            while((upper - lower) > sqrt(.Machine$double.eps) && 
                  (data_log_lik_prop < threshold)) {
                
                # increment the number of ESS steps
                increment_elem(tparam[[p]]$steps, k-1)
                
                # shrink the bracket
                if(theta < 0) {
                    lower <- theta
                } else {
                    upper <- theta
                }
                
                # sample a new point
                theta <- runif(1, lower, upper)
                
                # compute the proposal
                copy_vec(
                    dest = tparam[[p]]$draws_ess,
                    orig = cos(theta) * tparam[[p]]$draws_cur + 
                        sin(theta) * tparam[[p]]$draws_prop
                )
                
                # insert time-varying parameters
                insert_tparam(
                    tcovar = parmat,
                    values = 
                        tparam[[p]]$draws2par(
                            parameters = parmat[1,],
                            draws = tparam[[p]]$draws_ess),
                    col_ind = tparam[[p]]$col_ind,
                    tpar_inds = tparam[[p]]$tpar_inds_Cpp)    
                
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
                
            # if the bracket width is not equal to zero, update the draws, path, and dat log likelihood
            if((upper - lower) > sqrt(.Machine$double.eps)) {
                
                # transfer the draws
                copy_vec(dest = tparam[[p]]$draws_cur,
                         orig = tparam[[p]]$draws_ess)
                
                # copy the values
                copy_vec(dest = tparam[[p]]$tpar_cur,
                         orig = parmat[,tparam[[p]]$col_ind + 1])
                
                # copy the LNA path and the dat log likelihood
                copy_vec(dest = path$data_log_lik, orig = data_log_lik_prop)
                
                if(k == tparam_ess_control$n_updates & p == length(tparam)) {
                    copy_mat(dest = path$latent_path, orig = pathmat_prop)    
                }
                
                # record the final angle
                insert_elem(dest = tparam[[p]]$angles,
                            elem = theta,
                            ind  = k-1)
                
            } else {
                
                # recover the original time-varying parameters
                vec_2_mat(dest = parmat,
                          orig = tparam[[p]]$tpar_cur,
                          ind = tparam[[p]]$col_ind)
            }
        }
    }
        
    # update the ESS brackets
    if(iter != 0 && iter <= tparam_ess_control$bracket_update_iter) {
        
        for(p in seq_along(tparam)) {
            
            # angle residual
            copy_vec(
                dest = tparam[[p]]$angle_resid,
                orig = mean(tparam[[p]]$angles) - 
                    tparam[[p]]$angle_mean) 
            
            # angle variance
            copy_vec(
                dest = tparam[[p]]$angle_var,
                orig = tparam[[p]]$angle_resid^2 / (iter-1) + 
                    tparam[[p]]$angle_var * (iter-2) / (iter-1))
            
            # angle mean
            copy_vec(
                dest = tparam[[p]]$angle_mean,
                orig = mean(tparam[[p]]$angles) / (iter-1) + 
                    tparam[[p]]$angle_mean * (iter-2) / (iter-1))
            
            # set the new angle bracket
            if(iter == tparam_ess_control$bracket_update_iter) {
                copy_vec(
                    dest = tparam[[p]]$bracket_width,
                    orig = pmin(tparam_ess_control$bracket_scaling * 
                                    sqrt(tparam[[p]]$angle_var), 2*pi))    
            }
        }
    }
}
    
