#' Sample a new LNA path via elliptical slice sampling.
#'
#' @param path_cur list with the current LNA path along with its ODE paths
#' @param n_ess_updates number of elliptical slice sampling updates
#' @param svd_sqrt,svd_d,svd_U,svd_V objects for computing the SVD of LNA
#'   diffusion matrics
#' @param tparam_update if "joint" then time-varying parameters are updated jointly
#'   along with the LNA path
#' @inheritParams initialize_lna
#'
#' @return list with an updated LNA path along with its stochastic
#'   perturbations, the observed data log-likelihood, and the lna
#'   log-likelihood, and a record of the number of elliptical slice sampling
#'   proposals
#' @export
update_lna_path <-
        function(path_cur,
                 data,
                 lna_parameters,
                 lna_param_vec,
                 tparam,
                 pathmat_prop,
                 censusmat,
                 draws_prop,
                 ess_draws_prop,
                 emitmat,
                 flow_matrix,
                 stoich_matrix,
                 lna_times,
                 forcing_inds,
                 forcing_matrix,
                 lna_param_inds,
                 lna_const_inds,
                 lna_tcovar_inds,
                 lna_initdist_inds,
                 param_update_inds,
                 census_indices,
                 lna_event_inds,
                 measproc_indmat,
                 svd_sqrt,
                 svd_d,
                 svd_U,
                 svd_V,
                 lna_pointer,
                 lna_set_pars_pointer,
                 d_meas_pointer,
                 do_prevalence,
                 n_ess_updates,
                 tparam_update,
                 step_size) {
              
      copy_vec(path_cur$ess_record, rep(1, n_ess_updates))
      
      # perform the ESS updates, retaining the last state
      for(k in seq_len(n_ess_updates)) {

            # sample a new set of stochastic perturbations
            draw_normals2(ess_draws_prop)
            
            # choose a likelihood threshold
            threshold <- path_cur$data_log_lik + log(runif(1))
            
            # initial proposal, which also defines a bracket
            theta <- runif(1, 0, 2*pi)
            lower <- theta - 2*pi; upper <- theta
            
            # construct the first proposal
            copy_mat(draws_prop, cos(theta)*path_cur$draws + sin(theta)*ess_draws_prop)
            
            # propose time-varying parameter values if called for
            if(!is.null(tparam) && tparam_update == "joint") {
                  
                  # sample a new set of perturbations and construct the first proposal
                  for(p in seq_along(tparam)) {
                        
                        # sample perturbations
                        draw_normals(tparam[[p]]$draws_prop)
                        
                        # compute proposal
                        copy_vec(dest = tparam[[p]]$draws_ess, orig = cos(theta) * tparam[[p]]$draws_cur + sin(theta) * tparam[[p]]$draws_prop)
                        
                        # map to parameter
                        insert_tparam(tcovar    = lna_parameters,
                                      values    = tparam[[p]]$draws2par(parameters = lna_parameters[1,], draws = tparam[[p]]$draws_ess),
                                      col_ind   = tparam[[p]]$col_ind,
                                      tpar_inds = tparam[[p]]$tpar_inds)
                  }
            }
            
            # initialize the data log likelihood for the proposed path
            data_log_lik_prop <- NULL
            
            # map the perturbations to an LNA path
            try({
                    map_draws_2_lna(
                            pathmat           = pathmat_prop,
                            draws             = draws_prop,
                            lna_times         = lna_times,
                            lna_pars          = lna_parameters,
                            lna_param_vec     = lna_param_vec,
                            lna_param_inds    = lna_param_inds,
                            lna_tcovar_inds   = lna_tcovar_inds,
                            init_start        = lna_initdist_inds[1],
                            param_update_inds = param_update_inds,
                            stoich_matrix     = stoich_matrix,
                            forcing_inds      = forcing_inds,
                            forcing_matrix    = forcing_matrix,
                            svd_sqrt          = svd_sqrt,
                            svd_d             = svd_d,
                            svd_U             = svd_U,
                            svd_V             = svd_V,
                            lna_pointer       = lna_pointer,
                            set_pars_pointer  = lna_set_pars_pointer,
                            step_size         = step_size
                    )
            
                    census_lna(
                            path                = pathmat_prop,
                            census_path         = censusmat,
                            census_inds         = census_indices,
                            lna_event_inds      = lna_event_inds,
                            flow_matrix_lna     = flow_matrix,
                            do_prevalence       = do_prevalence,
                            init_state          = lna_parameters[1, lna_initdist_inds + 1, drop = TRUE],
                            forcing_matrix      = forcing_matrix
                    )
            
                    # evaluate the density of the incidence counts
                    evaluate_d_measure_LNA(
                            emitmat           = emitmat,
                            obsmat            = data,
                            censusmat         = censusmat,
                            measproc_indmat   = measproc_indmat,
                            lna_parameters    = lna_parameters,
                            lna_param_inds    = lna_param_inds,
                            lna_const_inds    = lna_const_inds,
                            lna_tcovar_inds   = lna_tcovar_inds,
                            param_update_inds = param_update_inds,
                            census_indices    = census_indices,
                            lna_param_vec     = lna_param_vec,
                            d_meas_ptr        = d_meas_pointer
                    )
            
                    # compute the data log likelihood
                    data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                    if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
            }, silent = TRUE)
            
            if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
            
            # continue proposing if not accepted
            while(!isTRUE(all.equal(lower, upper)) && (data_log_lik_prop < threshold)) {
            
                    # increment the number of ESS proposals for the current iteration
                    increment_elem(path_cur$ess_record, k-1)
                    
                    # shrink the bracket
                    if(theta < 0) {
                            lower <- theta
                    } else {
                            upper <- theta
                    }
            
                    # sample a new point
                    theta <- runif(1, lower, upper)
            
                    # construct the next LNA path proposal
                    copy_mat(draws_prop, cos(theta)*path_cur$draws + sin(theta)*ess_draws_prop)
                    
                    # construct the next proposal for time-varying parameters
                    if(!is.null(tparam) && tparam_update == "joint") {
                          
                          for(p in seq_along(tparam)) {
                                # compute proposal
                                copy_vec(dest = tparam[[p]]$draws_ess, orig = cos(theta) * tparam[[p]]$draws_cur + sin(theta) * tparam[[p]]$draws_prop)
                                
                                # map to parameter
                                insert_tparam(tcovar    = lna_parameters,
                                              values    = tparam[[p]]$draws2par(parameters = lna_parameters[1,], draws = tparam[[p]]$draws_ess),
                                              col_ind   = tparam[[p]]$col_ind,
                                              tpar_inds = tparam[[p]]$tpar_inds)
                          }
                    }
            
                    # initialize the data log likelihood for the proposed path
                    data_log_lik_prop <- NULL
            
                    # map the perturbations to an LNA path
                    try({
                          map_draws_2_lna(
                                pathmat           = pathmat_prop,
                                draws             = draws_prop,
                                lna_times         = lna_times,
                                lna_pars          = lna_parameters,
                                lna_param_vec     = lna_param_vec,
                                lna_param_inds    = lna_param_inds,
                                lna_tcovar_inds   = lna_tcovar_inds,
                                init_start        = lna_initdist_inds[1],
                                param_update_inds = param_update_inds,
                                stoich_matrix     = stoich_matrix,
                                forcing_inds      = forcing_inds,
                                forcing_matrix    = forcing_matrix,
                                svd_sqrt          = svd_sqrt,
                                svd_d             = svd_d,
                                svd_U             = svd_U,
                                svd_V             = svd_V,
                                lna_pointer       = lna_pointer,
                                set_pars_pointer  = lna_set_pars_pointer,
                                step_size         = step_size
                          )
            
                            census_lna(
                                    path                = pathmat_prop,
                                    census_path         = censusmat,
                                    census_inds         = census_indices,
                                    lna_event_inds      = lna_event_inds,
                                    flow_matrix_lna     = flow_matrix,
                                    do_prevalence       = do_prevalence,
                                    init_state          = lna_parameters[1, lna_initdist_inds + 1],
                                    forcing_matrix      = forcing_matrix
                            )
            
                            # evaluate the density of the incidence counts
                            evaluate_d_measure_LNA(
                                    emitmat           = emitmat,
                                    obsmat            = data,
                                    censusmat         = censusmat,
                                    measproc_indmat   = measproc_indmat,
                                    lna_parameters    = lna_parameters,
                                    lna_param_inds    = lna_param_inds,
                                    lna_const_inds    = lna_const_inds,
                                    lna_tcovar_inds   = lna_tcovar_inds,
                                    param_update_inds = param_update_inds,
                                    census_indices    = census_indices,
                                    lna_param_vec     = lna_param_vec,
                                    d_meas_ptr        = d_meas_pointer
                            )
            
                            # compute the data log likelihood
                            data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                            if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
                    }, silent = TRUE)
            
                    if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
            }
            
            # if the bracket width is not equal to zero, update the draws, path, and data log likelihood
            if(!isTRUE(all.equal(lower, upper))) {
            
                  # transfer the new path and residual path into the* sin(theta) path_prop list
                  copy_mat(path_cur$draws, draws_prop)
                  
                  if(!is.null(tparam) && tparam_update == "joint") {
                        # Copy the tparam draws
                        for(p in seq_along(tparam)) {
                              copy_vec(tparam[[p]]$draws_cur, tparam[[p]]$draws_ess)
                        }
                  }
            
                  # copy the LNA path and the data log likelihood
                  copy_mat(path_cur$lna_path, pathmat_prop)
                  copy_vec(path_cur$data_log_lik, data_log_lik_prop)
                    
            } else {
                  # if updating the time-varying parameters jointly with the LNA path
                  # recover the original time-varying parameter values
                  if(!is.null(tparam) && tparam_update == "joint") {
                        for(p in seq_along(tparam)) {
                              insert_tparam(tcovar    = lna_parameters,
                                            values    = tparam[[p]]$draws2par(parameters = lna_parameters[1,], draws = tparam[[p]]$draws_cur),
                                            col_ind   = tparam[[p]]$col_ind,
                                            tpar_inds = tparam[[p]]$tpar_inds)
                        }
                  }
            }
      }
}