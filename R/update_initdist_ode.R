#' Update initial compartment volumes via elliptical slice sampling.
#'
#' @param tparam list containing the time-varying parameters
#' @param path_cur list with the current path
#' @param n_ess_updates number of elliptical slice sampling updates
#' @inheritParams initialize_lna
#'
#' @return updated time-varying parameter values and lna path
#' @export
update_initdist_ode <-
        function(initdist_objects,
                 init_volumes_cur,
                 init_volumes_prop,
                 path_cur,
                 data,
                 ode_parameters,
                 ode_param_vec,
                 tparam,
                 pathmat_prop,
                 censusmat,
                 emitmat,
                 flow_matrix,
                 stoich_matrix,
                 ode_times,
                 forcing_inds,
                 forcing_matrix,
                 ode_param_inds,
                 ode_const_inds,
                 ode_tcovar_inds,
                 ode_initdist_inds,
                 param_update_inds,
                 census_indices,
                 ode_event_inds,
                 measproc_indmat,
                 ode_pointer,
                 ode_set_pars_pointer,
                 d_meas_pointer,
                 do_prevalence,
                 step_size,
                 initdist_ess,
                 initdist_bracket_width) {
              
      # initialize ess count
      ess_count <- 1
      
      # choose a likelihood threshold
      threshold <- path_cur$data_log_lik + log(runif(1))
      
      # initialize the data log likelihood for the proposed path
      data_log_lik_prop <- NULL
      
      # initial proposal, which also defines a bracket
      theta <- runif(1, 0, initdist_bracket_width)
      lower <- theta - initdist_bracket_width; upper <- theta
      
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
                           orig = cos(theta) * initdist_objects[[s]]$draws_cur + sin(theta) * initdist_objects[[s]]$draws_prop)
                  
                  # map to volumes
                  copy_vec2(dest = init_volumes_prop,
                            orig = initdist_objects[[s]]$comp_mean + 
                                  initdist_objects[[s]]$comp_sqrt_cov %*% initdist_objects[[s]]$draws_ess,
                            inds = initdist_objects[[s]]$comp_inds_Cpp)
                  
                  # check boundary conditions
                  bad_draws[s] <- 
                        any(init_volumes_prop[initdist_objects[[s]]$comp_inds_R] < 0) | 
                        any(init_volumes_prop[initdist_objects[[s]]$comp_inds_R] > initdist_objects[[s]]$comp_size)
            }
      }
      
      # check boundary conditions, if valid integrate odes
      if(any(bad_draws)) {
            data_log_lik_prop <- -Inf 
            
      } else {
            
            # copy the proposed state into the parameter matrix
            pars2lnapars2(ode_parameters, init_volumes_prop, ode_initdist_inds[1])
            
            # map time varying parameters if called for
            if(!is.null(tparam)) {
                  for(p in seq_along(tparam)) {
                        # map to parameter
                        insert_tparam(tcovar    = ode_parameters,
                                      values    = tparam[[p]]$draws2par(parameters = ode_parameters[1,], draws = tparam[[p]]$draws_cur),
                                      col_ind   = tparam[[p]]$col_ind,
                                      tpar_inds = tparam[[p]]$tpar_inds)
                  }      
            }
            
            # map the perturbations to an LNA path
            try({
                  map_pars_2_ode(
                        pathmat           = pathmat_prop,
                        ode_times         = ode_times,
                        ode_pars          = ode_parameters,
                        init_start        = ode_initdist_inds[1],
                        ode_param_inds    = ode_param_inds,
                        ode_tcovar_inds   = ode_tcovar_inds,
                        param_update_inds = param_update_inds,
                        stoich_matrix     = stoich_matrix,
                        forcing_inds      = forcing_inds,
                        forcing_matrix    = forcing_matrix,
                        ode_pointer       = ode_pointer,
                        set_pars_pointer  = ode_set_pars_pointer,
                        step_size         = step_size
                  )
                  
                  census_lna(
                        path                = pathmat_prop,
                        census_path         = censusmat,
                        census_inds         = census_indices,
                        lna_event_inds      = ode_event_inds,
                        flow_matrix_lna     = flow_matrix,
                        do_prevalence       = do_prevalence,
                        init_state          = init_volumes_prop,
                        forcing_matrix      = forcing_matrix
                  )
                  
                  # evaluate the density of the incidence counts
                  evaluate_d_measure_LNA(
                        emitmat           = emitmat,
                        obsmat            = data,
                        censusmat         = censusmat,
                        measproc_indmat   = measproc_indmat,
                        lna_parameters    = ode_parameters,
                        lna_param_inds    = ode_param_inds,
                        lna_const_inds    = ode_const_inds,
                        lna_tcovar_inds   = ode_tcovar_inds,
                        param_update_inds = param_update_inds,
                        census_indices    = census_indices,
                        lna_param_vec     = ode_param_vec,
                        d_meas_ptr        = d_meas_pointer
                  )
                  
                  # compute the data log likelihood
                  data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                  if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
            }, silent = TRUE)
            
            if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
      }

      # continue proposing if not accepted
      while((upper - lower) > sqrt(.Machine$double.eps) && (data_log_lik_prop < threshold)) {
      
              # increment the number of ESS proposals for the current iteration
              ess_count <- ess_count + 1
      
              # shrink the bracket
              if(theta < 0) {
                      lower <- theta
              } else {
                      upper <- theta
              }
      
              # sample a new point
              theta <- runif(1, lower, upper)
              
              # construct the next proposal
              for(s in seq_along(initdist_objects)) {
                    
                    # if the state is not fixed draw new values
                    if(!initdist_objects[[s]]$fixed) {
                          
                          # compute the linear combination
                          copy_vec(dest = initdist_objects[[s]]$draws_ess, 
                                   orig = cos(theta) * initdist_objects[[s]]$draws_cur + sin(theta) * initdist_objects[[s]]$draws_prop)
                          
                          # map to volumes
                          copy_vec2(dest = init_volumes_prop,
                                    orig = initdist_objects[[s]]$comp_mean + 
                                          initdist_objects[[s]]$comp_sqrt_cov %*% initdist_objects[[s]]$draws_ess,
                                    inds = initdist_objects[[s]]$comp_inds_Cpp)
                          
                          # check boundary conditions
                          bad_draws[s] <- 
                                any(init_volumes_prop[initdist_objects[[s]]$comp_inds_R] < 0) | 
                                any(init_volumes_prop[initdist_objects[[s]]$comp_inds_R] > initdist_objects[[s]]$comp_size)
                    }
              }
      
              # check boundary conditions, if valid integrate odes
              if(any(bad_draws)) {
                    data_log_lik_prop <- -Inf 
                    
              } else {
                    
                    # copy the proposed state into the parameter matrix
                    pars2lnapars2(ode_parameters, init_volumes_prop, ode_initdist_inds[1])
                    
                    # map time varying parameters if called for
                    if(!is.null(tparam)) {
                          for(p in seq_along(tparam)) {
                                # map to parameter
                                insert_tparam(tcovar    = ode_parameters,
                                              values    = tparam[[p]]$draws2par(parameters = ode_parameters[1,], draws = tparam[[p]]$draws_cur),
                                              col_ind   = tparam[[p]]$col_ind,
                                              tpar_inds = tparam[[p]]$tpar_inds)
                          }      
                    }
                    
                    # map the perturbations to an LNA path
                    try({
                          map_pars_2_ode(
                                pathmat           = pathmat_prop,
                                ode_times         = ode_times,
                                ode_pars          = ode_parameters,
                                init_start        = ode_initdist_inds[1],
                                ode_param_inds    = ode_param_inds,
                                ode_tcovar_inds   = ode_tcovar_inds,
                                param_update_inds = param_update_inds,
                                stoich_matrix     = stoich_matrix,
                                forcing_inds      = forcing_inds,
                                forcing_matrix    = forcing_matrix,
                                ode_pointer       = ode_pointer,
                                set_pars_pointer  = ode_set_pars_pointer,
                                step_size         = step_size
                          )
                          
                          census_lna(
                                path                = pathmat_prop,
                                census_path         = censusmat,
                                census_inds         = census_indices,
                                lna_event_inds      = ode_event_inds,
                                flow_matrix_lna     = flow_matrix,
                                do_prevalence       = do_prevalence,
                                init_state          = init_volumes_prop,
                                forcing_matrix      = forcing_matrix
                          )
                          
                          # evaluate the density of the incidence counts
                          evaluate_d_measure_LNA(
                                emitmat           = emitmat,
                                obsmat            = data,
                                censusmat         = censusmat,
                                measproc_indmat   = measproc_indmat,
                                lna_parameters    = ode_parameters,
                                lna_param_inds    = ode_param_inds,
                                lna_const_inds    = ode_const_inds,
                                lna_tcovar_inds   = ode_tcovar_inds,
                                param_update_inds = param_update_inds,
                                census_indices    = census_indices,
                                lna_param_vec     = ode_param_vec,
                                d_meas_ptr        = d_meas_pointer
                          )
                          
                          # compute the data log likelihood
                          data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                          if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
                    }, silent = TRUE)
                    
                    if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
              }
      }
      
      # if the bracket width is not equal to zero, update the draws, path, and data log likelihood
      if((upper - lower) > sqrt(.Machine$double.eps)) {
      
              # transfer the new initial volumes and draws
              for(s in seq_along(initdist_objects)) {
                    if(!initdist_objects[[s]]$fixed) {
                          
                          # copy the N(0,1) draws
                          copy_vec(dest = initdist_objects[[s]]$draws_cur,
                                   orig = initdist_objects[[s]]$draws_ess)
                          
                          # copy the initial compartment volumes
                          copy_vec2(dest = init_volumes_cur,
                                    orig = init_volumes_prop[initdist_objects[[s]]$comp_inds_R],
                                    inds = initdist_objects[[s]]$comp_inds_Cpp)
                    }
              }
      
              # copy the ode path and the data log likelihood
              copy_mat(path_cur$ode_path, pathmat_prop)
              copy_vec(path_cur$data_log_lik, data_log_lik_prop)
              
      } else {
            
            # insert the original initial compartment counts into the parameter matrix
            pars2lnapars2(ode_parameters, init_volumes_cur, ode_initdist_inds[1])
            
            # recover the original time-varying parameter values
            if(!is.null(tparam)) {
                  for(p in seq_along(tparam)) {
                        insert_tparam(tcovar    = ode_parameters,
                                      values    = tparam[[p]]$draws2par(parameters = ode_parameters[1,], draws = tparam[[p]]$draws_cur),
                                      col_ind   = tparam[[p]]$col_ind,
                                      tpar_inds = tparam[[p]]$tpar_inds)
                  }     
            }
      }
      
      copy_vec(initdist_ess, ess_count)
}