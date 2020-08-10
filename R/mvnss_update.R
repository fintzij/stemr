#' Update model parameters via factor slice sampling
#'
#' @inheritParams mvnmh_update
#'
#' @return update the model parameters, path, and likelihood
#' @export
mvnss_update <- 
    function(param_blocks, 
             ind,
             iter,
             parmat, 
             dat,
             path,
             pathmat_prop,
             tparam,
             times,
             flow_matrix,
             stoich_matrix,
             censusmat,
             emitmat,
             param_vec,
             initdist_inds,
             census_indices,
             event_inds,
             measproc_indmat,
             forcing_inds,
             forcing_tcov_inds,
             forcings_out,
             forcing_transfers,
             proc_pointer,
             d_meas_pointer,
             do_prevalence,
             step_size,
             svd_d = NULL,
             svd_U = NULL,
             svd_V = NULL) {
        
        # sample the likelihood threshold
        threshold <- path$data_log_lik + param_blocks[[ind]]$log_pd - rexp(1)
        
        # sample the hit-and-run direction
        if(param_blocks[[ind]]$nugget_sequence[iter] != 0) {
            
            draw_normals(param_blocks[[ind]]$mvnss_objects$mvn_direction)
            sample_unit_sphere(param_blocks[[ind]]$mvnss_objects$har_direction)
            
            # compute the proposal
            copy_vec(dest = param_blocks[[ind]]$mvnss_objects$mvnss_propvec, 
                     orig = normalise2((1 - param_blocks[[ind]]$nugget_sequence[iter]) * 
                                           normalise2(param_blocks[[ind]]$mvnss_objects$mvn_direction %*% 
                                                          param_blocks[[ind]]$kernel_cov_chol, 2) + 
                                           param_blocks[[ind]]$nugget_sequence[iter] * 
                                           param_blocks[[ind]]$mvnss_objects$har_direction, 2))
            
        } else {
            draw_normals(param_blocks[[ind]]$mvnss_objects$mvn_direction)
            copy_vec(dest = param_blocks[[ind]]$mvnss_objects$mvnss_propvec, 
                     orig = normalise2(param_blocks[[ind]]$mvnss_objects$mvn_direction %*% 
                                           param_blocks[[ind]]$kernel_cov_chol, 2))
        }
        
        # construct the approximate bracket
        center <- runif(1)
        lower  <- -param_blocks[[ind]]$mvnss_objects$bracket_width * center
        upper  <- lower + param_blocks[[ind]]$mvnss_objects$bracket_width
        
        # initialize the log-posterior at the endpoints
        logpost_lower <- NULL
        logpost_upper <- NULL
        logpost_prop  <- -Inf 
        
        # # step out lower bound
        while(is.null(logpost_lower) || threshold < logpost_lower) {
            
            # lower end of the bracket on the estimation scale
            copy_vec(dest = param_blocks[[ind]]$pars_prop_est, 
                     orig = param_blocks[[ind]]$pars_est + 
                         lower * param_blocks[[ind]]$mvnss_objects$mvnss_propvec)
            
            # get the parameters on the natural scale
            copy_vec(param_blocks[[ind]]$pars_prop_nat, 
                     param_blocks[[ind]]$priors$from_estimation_scale(
                         param_blocks[[ind]]$pars_prop_est))
            
            # compute the prior density
            logprior_lower <- 
                param_blocks[[ind]]$priors$logprior(param_blocks[[ind]]$pars_prop_est)
            
            # if the log prior is not -Inf, find the path
            if(logprior_lower != -Inf) {
                
                # insert parameters into the parameter proposal matrix
                pars2parmat(parmat  = parmat,
                            pars    = param_blocks[[ind]]$pars_prop_nat,
                            colinds = param_blocks[[ind]]$param_inds_Cpp)
                
                # compute the time-varying parameters if necessary
                if(!is.null(tparam)) {
                    for(p in seq_along(tparam)) {
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
                
                # initialize data log likelihood
                loglik_lower <- NULL
                
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
                    loglik_lower <- sum(emitmat[, -1][measproc_indmat])
                    
                    if (is.nan(loglik_lower)) loglik_lower <- -Inf
                }, silent = TRUE)
                
                if(is.null(loglik_lower)) loglik_lower <- -Inf      
                
            } else {
                loglik_lower <- -Inf
            }
            
            # compute log-posterior
            logpost_lower <- loglik_lower + logprior_lower
            
            # step out the bracket if necessary
            if(threshold < logpost_lower) {
                
                # decrease the lower endpoint of the bracket
                lower <- lower - param_blocks[[ind]]$mvnss_objects$bracket_width
                
                # increment the number of expansions
                increment_elem(param_blocks[[ind]]$mvnss_objects$n_expansions, 0)
            }
        }
        
        # step out upper
        while(is.null(logpost_upper) || threshold < logpost_upper) {
            
            # upper end of the bracket on the estimation scale
            copy_vec(dest = param_blocks[[ind]]$pars_prop_est, 
                     orig = param_blocks[[ind]]$pars_est + 
                         upper * param_blocks[[ind]]$mvnss_objects$mvnss_propvec)
            
            # get the parameters on the natural scale
            copy_vec(param_blocks[[ind]]$pars_prop_nat, 
                     param_blocks[[ind]]$priors$from_estimation_scale(
                         param_blocks[[ind]]$pars_prop_est))
            
            # compute the prior density
            logprior_upper <-
                param_blocks[[ind]]$priors$logprior(param_blocks[[ind]]$pars_prop_est)
            
            # if the log prior is not -Inf, find the path
            if(logprior_upper != -Inf) {
                
                # insert parameters into the parameter proposal matrix
                pars2parmat(parmat  = parmat,
                            pars    = param_blocks[[ind]]$pars_prop_nat,
                            colinds = param_blocks[[ind]]$param_inds_Cpp)
                
                # compute the time-varying parameters if necessary
                if(!is.null(tparam)) {
                    for(p in seq_along(tparam)) {
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
                
                # initialize data log likelihood
                loglik_upper <- NULL
                
                # map the perturbations to an LNA path
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
                    loglik_upper <- sum(emitmat[, -1][measproc_indmat])
                    
                    if (is.nan(loglik_upper)) loglik_upper <- -Inf
                }, silent = TRUE)
                
                if(is.null(loglik_upper)) loglik_upper <- -Inf      
                
            } else {
                loglik_upper <- -Inf
            }
            
            # compute log-posterior
            logpost_upper <- loglik_upper + logprior_upper
            
            # step out the bracket if necessary
            if(threshold < logpost_upper) {
                
                # decrease the upper endpoint of the bracket
                upper <- upper + param_blocks[[ind]]$mvnss_objects$bracket_width
                
                # increment the number of expansions
                increment_elem(param_blocks[[ind]]$mvnss_objects$n_expansions, 0)
            }
        }
        
        # sample from the bracket
        while((upper - lower) > sqrt(.Machine$double.eps) && (logpost_prop < threshold)) {
            
            # sample uniformly in the bracket
            prop <- runif(1, lower, upper)
            
            # proposal on the estimation scale
            copy_vec(dest = param_blocks[[ind]]$pars_prop_est, 
                     orig = param_blocks[[ind]]$pars_est + 
                         prop * param_blocks[[ind]]$mvnss_objects$mvnss_propvec)
            
            # get the parameters on the natural scale
            copy_vec(param_blocks[[ind]]$pars_prop_nat, 
                     param_blocks[[ind]]$priors$from_estimation_scale(
                         param_blocks[[ind]]$pars_prop_est))
            
            # compute the prior density
            logprior_prop <-
                param_blocks[[ind]]$priors$logprior(param_blocks[[ind]]$pars_prop_est)
            
            # if the log prior is not -Inf, find the path
            if(logprior_prop != -Inf) {
                
                # insert parameters into the parameter proposal matrix
                pars2parmat(parmat  = parmat,
                            pars    = param_blocks[[ind]]$pars_prop_nat,
                            colinds = param_blocks[[ind]]$param_inds_Cpp)
                
                # compute the time-varying parameters if necessary
                if(!is.null(tparam)) {
                    for(p in seq_along(tparam)) {
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
                
                # initialize data log likelihood
                loglik_prop <- NULL
                
                # map the perturbations to an LNA path
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
                    loglik_prop <- sum(emitmat[, -1][measproc_indmat])
                    
                    if (is.nan(loglik_prop)) loglik_prop <- -Inf
                }, silent = TRUE)
                
                if(is.null(loglik_prop)) loglik_prop <- -Inf      
                
            } else {
                loglik_prop <- -Inf
            }
            
            # compute log-posterior
            logpost_prop <- loglik_prop + logprior_prop
            
            # shrink the bracket if necessary
            if(threshold > logpost_prop) {
                
                # adjust the bounds
                if(prop < 0) {
                    lower <- prop
                } else {
                    upper <- prop
                }
                
                # increment the number of contractions 
                increment_elem(param_blocks[[ind]]$mvnss_objects$n_contractions, 0)
            }
        }
        
        if((upper - lower) > sqrt(.Machine$double.eps)) {
            
            # update log likelihood and prior
            copy_vec(path$data_log_lik, loglik_prop)  
            copy_vec(param_blocks[[ind]]$log_pd, logprior_prop) 
            
            # copy parameters
            copy_vec(param_blocks[[ind]]$pars_nat, param_blocks[[ind]]$pars_prop_nat)
            copy_vec(param_blocks[[ind]]$pars_est, param_blocks[[ind]]$pars_prop_est)
            
            # copy time-varying parameters
            if(!is.null(tparam)) {
                for(p in seq_along(tparam)) {
                        copy_vec(dest = tparam[[p]]$tpar_cur,
                                 orig = parmat[,tparam[[p]]$col_ind + 1])    
                }
            }
            
            # copy the path matrix
            copy_mat(path$latent_path, pathmat_prop)
            
        } else {
            
            # need to reset the parmat matrix
            pars2parmat(parmat  = parmat,
                        pars    = param_blocks[[ind]]$pars_nat,
                        colinds = param_blocks[[ind]]$param_inds_Cpp)
            
            if(!is.null(tparam)) {
                for(p in seq_along(tparam)) {
                    vec_2_mat(dest = parmat,
                              orig = tparam[[p]]$tpar_cur,
                              ind = tparam[[p]]$col_ind)
                }
            }
        }
        
        # adapt the MCMC kernel
        if(iter < param_blocks[[ind]]$control$stop_adaptation) {
            
            # calculate the residual
            copy_vec(param_blocks[[ind]]$kernel_resid,
                     param_blocks[[ind]]$pars_est - param_blocks[[ind]]$kernel_mean)
            
            # update the empirical covariance matrix
            copy_mat(param_blocks[[ind]]$kernel_cov,
                     param_blocks[[ind]]$kernel_cov + 
                         param_blocks[[ind]]$gain_factors[iter] * 
                         (param_blocks[[ind]]$kernel_resid %o% param_blocks[[ind]]$kernel_resid - 
                              param_blocks[[ind]]$kernel_cov))
            
            # update the empirical mean
            copy_vec(param_blocks[[ind]]$kernel_mean,
                     param_blocks[[ind]]$kernel_mean + 
                         param_blocks[[ind]]$gain_factors[iter] * 
                         param_blocks[[ind]]$kernel_resid)
            
            # compute the cholesky
            comp_chol(param_blocks[[ind]]$kernel_cov_chol, 
                      param_blocks[[ind]]$kernel_cov)
        }
        
        # adapt the bracket width
        param_blocks[[ind]]$mvnss_objects$bracket_width = 
            max(param_blocks[[ind]]$control$bracket_limits[1],
                min(param_blocks[[ind]]$control$bracket_limits[2],
                    exp(log(param_blocks[[ind]]$mvnss_objects$bracket_width) +
                            param_blocks[[ind]]$mvnss_objects$bracket_adaptations[iter] * 
                            (param_blocks[[ind]]$mvnss_objects$n_expansions / 
                                 (param_blocks[[ind]]$mvnss_objects$n_expansions + 
                                      param_blocks[[ind]]$mvnss_objects$n_contractions) - 0.5))))
    }
