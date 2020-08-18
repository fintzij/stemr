#' Multivariate normal Metropolis-Hastings update
#'
#' @param param_blocks list of parameter blocks
#' @param ind index of parameter block to be updated
#' @param iter MCMC iteration
#' @param params_cur matrix with current parameters
#' @param params_prop matrix with proposed parameters
#' @param dat data matrix
#' @param path path list
#' @param pathmat_prop matrix for proposed paths
#' @param tparam list of time-varying parameters
#' @param times vector of census times
#' @param flow_matrix flow matrix
#' @param stoich_matrix stoichiometry matrix
#' @param censusmat census matrix
#' @param emitmat emission matrix
#' @param param_vec parameter vector for use in emission distribution
#' @param initdist_inds indices of initial compartment counts
#' @param census_indices indices where the path should be censused
#' @param event_inds event indices
#' @param measproc_indmat measurement process indices
#' @param forcing_inds forcing indices
#' @param forcing_tcov_inds time varying covariate forcings
#' @param forcings_out matrix with outflow from forcings
#' @param forcing_transfers transfer matrix 
#' @param proc_pointer C++ pointer for latent process
#' @param d_meas_pointer C++ pointer for emission distribution
#' @param do_prevalence should prevalence be computed
#' @param step_size initial step size for ODE solvers
#' @param svd_d,svd_U,svd_V SVD objects for LNA, NULL if using the ODE approx
#'
#' @return update the model parameters, path, and likelihood
#' @export
mvnmh_update = 
    function(param_blocks, 
             ind,
             iter,
             parmat,
             dat,
             path,
             pathmat_prop,
             tparam,
             census_times,
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
        
        # propose new parameter values
        propose_mvnmh(
            params_prop = param_blocks[[ind]]$pars_prop_est,
            params_cur = param_blocks[[ind]]$pars_est,
            kernel_cov_chol = param_blocks[[ind]]$kernel_cov_chol,
            nugget = param_blocks[[ind]]$nugget_sequence[iter]
        )
        
        # go back to the natural scale
        copy_vec(param_blocks[[ind]]$pars_prop_nat, 
                 param_blocks[[ind]]$priors$from_estimation_scale(
                     param_blocks[[ind]]$pars_prop_est)
        )
        
        # compute the log prior for the proposed parameters
        param_blocks[[ind]]$log_pd_prop = 
            param_blocks[[ind]]$priors$logprior(
                param_blocks[[ind]]$pars_prop_est
        )
        
        # insert parameters into the parameter proposal matrix
        pars2parmat(parmat = parmat,
                    pars = param_blocks[[ind]]$pars_prop_nat,
                    colinds = param_blocks[[ind]]$param_inds_Cpp)
        
        # update time-varying parameters if necessary
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
        
        # set the data log likelihood for the proposal to NULL
        data_log_lik_prop <- NULL
            
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
        
        # if integration failed then reject    
        if (is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
        
        ## Compute the acceptance probability
        acceptance_prob <- 
            (data_log_lik_prop + param_blocks[[ind]]$log_pd_prop) - 
            (path$data_log_lik + param_blocks[[ind]]$log_pd)
        
        # Accept/Reject via metropolis-hastings
        if (acceptance_prob >= min(0, log(runif(1)))) {
            
            ### ACCEPTANCE
            increment_elem(param_blocks[[ind]]$mvnmh_objects$acceptances, 0)
            
            # update log likelihood and prior
            copy_vec(dest = path$data_log_lik, 
                     orig = data_log_lik_prop)  
            copy_vec(dest = param_blocks[[ind]]$log_pd, 
                     orig = param_blocks[[ind]]$log_pd_prop) 
            
            # copy parameters
            copy_vec(dest = param_blocks[[ind]]$pars_nat, 
                     orig = param_blocks[[ind]]$pars_prop_nat)
            copy_vec(dest = param_blocks[[ind]]$pars_est, 
                     orig = param_blocks[[ind]]$pars_prop_est)
            
            # copy time-varying parameters
            if(!is.null(tparam)) {
                for(p in seq_along(tparam)) {
                    if(tparam[[p]]$init_dep) {
                        copy_vec(dest = tparam[[p]]$tpar_cur,
                                 orig = parmat[,tparam[[p]]$col_ind + 1])    
                    }
                }
            }
            
            # copy latent path
            copy_pathmat(path$latent_path, pathmat_prop)
            
        } else {
            
            # need to reset the params_prop matrix
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
        if (iter <= param_blocks[[ind]]$control$stop_adaptation) {
            
            # Adapt the proposal kernel
            copy_vec(dest = param_blocks[[ind]]$mvnmh_objects$proposal_scaling,
                     orig = min(exp(log(param_blocks[[ind]]$mvnmh_objects$proposal_scaling) +
                                 param_blocks[[ind]]$gain_factors[iter] * 
                                 (min(exp(acceptance_prob), 1) - 
                                      param_blocks[[ind]]$control$target_acceptance)),
                         param_blocks[[ind]]$control$max_scaling))
            
            # calculate the residual
            copy_vec(dest = param_blocks[[ind]]$kernel_resid,
                     orig = param_blocks[[ind]]$pars_est - param_blocks[[ind]]$kernel_mean)
            
            # update the empirical covariance matrix
            copy_mat(dest = param_blocks[[ind]]$kernel_cov,
                     orig = param_blocks[[ind]]$kernel_cov + 
                         param_blocks[[ind]]$gain_factors[iter] * 
                         (param_blocks[[ind]]$kernel_resid %o% param_blocks[[ind]]$kernel_resid - 
                              param_blocks[[ind]]$kernel_cov))
            
            # update the empirical mean
            copy_vec(dest = param_blocks[[ind]]$kernel_mean,
                     orig = param_blocks[[ind]]$kernel_mean + 
                         param_blocks[[ind]]$gain_factors[iter] * 
                         param_blocks[[ind]]$kernel_resid)
            
            # compute the cholesky
            comp_chol(C = param_blocks[[ind]]$kernel_cov_chol, 
                      M = param_blocks[[ind]]$mvnmh_objects$proposal_scaling * 
                          param_blocks[[ind]]$kernel_cov)
        }
    }
