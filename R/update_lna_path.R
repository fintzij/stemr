#' Sample a new LNA path via elliptical slice sampling.
#'
#' @param path_cur list with the current LNA path along with its ODE paths
#' @param n_ess_updates number of elliptical slice sampling updates
#' @param svd_d,svd_U,svd_V objects for computing the SVD of LNA
#'   diffusion matrics
#' @param tparam_update if TRUE then time-varying parameters are updated jointly
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
                 init_volumes_cur,
                 init_volumes_prop,
                 initdist_objects,
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
                 forcing_tcov_inds,
                 forcings_out,
                 forcing_transfers,
                 lna_param_inds,
                 lna_const_inds,
                 lna_tcovar_inds,
                 lna_initdist_inds,
                 param_update_inds,
                 census_indices,
                 lna_event_inds,
                 measproc_indmat,
                 svd_d,
                 svd_U,
                 svd_V,
                 lna_pointer,
                 lna_set_pars_pointer,
                 d_meas_pointer,
                 do_prevalence,
                 n_ess_updates,
                 ess_schedule,
                 lna_bracket_width,
                 joint_tparam_update,
                 joint_initdist_update,
                 joint_strata_update,
                 step_size) {
              
      step_count <- matrix(1.0, nrow = n_ess_updates, ncol = length(ess_schedule[[1]]))
      ess_angles <- matrix(1.0, nrow = n_ess_updates, ncol = length(ess_schedule[[1]]))

      # perform the ESS updates, retaining the last state
      for(k in seq_len(n_ess_updates)) {
            
            # sample the order of the updates
            ess_order <- sample.int(length(ess_schedule[[1]]), replace = FALSE)
            
            # do those updates!
            for(j in ess_order) {
                  
                  # sample a new set of stochastic perturbations
                  ess_draws_prop[ess_schedule[[1]][[j]],] <- rnorm(ess_draws_prop[ess_schedule[[1]][[j]],])
                  
                  # choose a likelihood threshold
                  threshold <- path_cur$data_log_lik + log(runif(1))
                  
                  # initial proposal, which also defines a bracket
                  # theta <- runif(1, 0, ess_bracket_width)
                  # lower <- theta - ess_bracket_width; upper <- theta
                  pos <- runif(1) 
                  lower <- -lna_bracket_width * pos; upper <- lower + lna_bracket_width
                  theta <- runif(1, lower, upper)
                  
                  # initialize the data log likelihood for the proposed path
                  data_log_lik_prop <- NULL
                  
                  # propose a new initial state
                  if(joint_initdist_update) {
                        
                        bad_draws <- vector("logical", length(initdist_objects))
                        
                        #################################################################
                        ### Resume here - update initdists on same schedule as strata ###
                        #################################################################
                        
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
                  }
                  
                  if(joint_initdist_update && any(bad_draws)) {
                        data_log_lik_prop <- -Inf
                        
                  } else {
                        
                        if(joint_initdist_update) {
                              # copy the initial volumes into the parameter matrix
                              pars2lnapars2(lna_parameters, init_volumes_prop, lna_initdist_inds[1])
                        }
                        
                        # construct the first proposal
                        copy_mat(dest = draws_prop, orig = cos(theta) * path_cur$draws + sin(theta) * ess_draws_prop)
                        
                        # propose time-varying parameter values if called for
                        if(!is.null(tparam)) {
                              
                              if(joint_tparam_update) {
                                    
                                    # sample a new set of perturbations and construct the first proposal
                                    for(p in seq_along(tparam)) {
                                          
                                          # sample perturbations
                                          draw_normals(tparam[[p]]$draws_prop)
                                          
                                          # compute proposal
                                          copy_vec(dest = tparam[[p]]$draws_ess, 
                                                   orig = cos(theta) * tparam[[p]]$draws_cur + sin(theta) * tparam[[p]]$draws_prop)
                                          
                                          # map to parameter
                                          insert_tparam(tcovar    = lna_parameters,
                                                        values    = tparam[[p]]$draws2par(parameters = lna_parameters[1,], 
                                                                                          draws = tparam[[p]]$draws_ess),
                                                        col_ind   = tparam[[p]]$col_ind,
                                                        tpar_inds = tparam[[p]]$tpar_inds)
                                    }
                                    
                              } else {
                                    # or recompute time varying parameters 
                                    
                                    for(p in seq_along(tparam)) {
                                          # map to parameter
                                          insert_tparam(tcovar    = lna_parameters,
                                                        values    = tparam[[p]]$draws2par(parameters = lna_parameters[1,],
                                                                                          draws = tparam[[p]]$draws_cur),
                                                        col_ind   = tparam[[p]]$col_ind,
                                                        tpar_inds = tparam[[p]]$tpar_inds)
                                    }
                              }
                        }
                        
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
                                    forcing_tcov_inds = forcing_tcov_inds,
                                    forcings_out      = forcings_out,
                                    forcing_transfers = forcing_transfers,
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
                                    lna_pars            = lna_parameters,
                                    forcing_inds        = forcing_inds,
                                    forcing_tcov_inds   = forcing_tcov_inds,
                                    forcings_out        = forcings_out,
                                    forcing_transfers   = forcing_transfers
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
                  
                  # continue proposing if not accepted
                  while((upper - lower) > sqrt(.Machine$double.eps) && (data_log_lik_prop < threshold)) {
                        
                        # increment the number of ESS proposals for the current iteration
                        step_count[k] <- step_count[k] + 1
                        
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
                        }
                        
                        if(joint_initdist_update && any(bad_draws)) {
                              data_log_lik_prop <- -Inf
                              
                        } else {
                              
                              # copy the initial volumes into the parameter matrix
                              if(joint_initdist_update) {
                                    pars2lnapars2(lna_parameters, init_volumes_prop, lna_initdist_inds[1])
                              }
                              
                              # construct the next path proposal
                              copy_mat(dest = draws_prop, orig = cos(theta) * path_cur$draws + sin(theta) * ess_draws_prop)
                              
                              # propose time-varying parameter values if called for
                              if(!is.null(tparam)) {
                                    
                                    if(joint_tparam_update) {
                                          
                                          # sample a new set of perturbations and construct the first proposal
                                          for(p in seq_along(tparam)) {
                                                
                                                # compute proposal
                                                copy_vec(dest = tparam[[p]]$draws_ess, 
                                                         orig = cos(theta) * tparam[[p]]$draws_cur + sin(theta) * tparam[[p]]$draws_prop)
                                                
                                                # map to parameter
                                                insert_tparam(tcovar    = lna_parameters,
                                                              values    = tparam[[p]]$draws2par(parameters = lna_parameters[1,], 
                                                                                                draws = tparam[[p]]$draws_ess),
                                                              col_ind   = tparam[[p]]$col_ind,
                                                              tpar_inds = tparam[[p]]$tpar_inds)
                                          }
                                          
                                    } else {
                                          # or recompute time varying parameters 
                                          
                                          for(p in seq_along(tparam)) {
                                                # map to parameter
                                                insert_tparam(tcovar    = lna_parameters,
                                                              values    = tparam[[p]]$draws2par(parameters = lna_parameters[1,],
                                                                                                draws = tparam[[p]]$draws_cur),
                                                              col_ind   = tparam[[p]]$col_ind,
                                                              tpar_inds = tparam[[p]]$tpar_inds)
                                          }
                                    }
                              }
                              
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
                                          forcing_tcov_inds = forcing_tcov_inds,
                                          forcings_out      = forcings_out,
                                          forcing_transfers = forcing_transfers,
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
                                          lna_pars            = lna_parameters,
                                          forcing_inds        = forcing_inds,
                                          forcing_tcov_inds   = forcing_tcov_inds,
                                          forcings_out        = forcings_out,
                                          forcing_transfers   = forcing_transfers
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
                  }
                  
                  # if the bracket width is not equal to zero, update the draws, path, and data log likelihood
                  if((upper - lower) > sqrt(.Machine$double.eps)) {
                        
                        # transfer the new initial volumes and draws (volumes already in parameter matrix)
                        if(joint_initdist_update) {
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
                        }
                        
                        # copy time varying parameter draws (mapped values already in parameter matrix)
                        if(!is.null(tparam) && joint_tparam_update) {
                              for(p in seq_along(tparam)) {
                                    copy_vec(dest = tparam[[p]]$draws_cur, orig = tparam[[p]]$draws_ess)
                              }
                        }
                        
                        # transfer the new path and residual path into the* sin(theta) path_prop list
                        copy_mat(dest = path_cur$draws, orig = draws_prop)
                        
                        # copy the LNA path and the data log likelihood
                        copy_mat(dest = path_cur$lna_path, orig = pathmat_prop)
                        copy_vec(dest = path_cur$data_log_lik, orig = data_log_lik_prop)
                        
                        # record the final angle
                        ess_angles[k] <- theta
                        
                  } else {
                        
                        # insert the original compartment counts back into the parameter matrix
                        if(joint_initdist_update) {
                              pars2lnapars2(lna_parameters, init_volumes_cur, lna_initdist_inds[1])
                        }
                        
                        # if updating the time-varying parameters jointly with the LNA path
                        # recover the original time-varying parameter values
                        if(!is.null(tparam) && (joint_tparam_update || joint_initdist_update)) {
                              for(p in seq_along(tparam)) {
                                    insert_tparam(tcovar    = lna_parameters,
                                                  values    = tparam[[p]]$draws2par(parameters = lna_parameters[1,], 
                                                                                    draws = tparam[[p]]$draws_cur),
                                                  col_ind   = tparam[[p]]$col_ind,
                                                  tpar_inds = tparam[[p]]$tpar_inds)
                              }
                        }
                  }      
            }
      }
      
      # copy the step and angle records
      copy_vec(dest = path_cur$step_record, orig = step_count)
      copy_vec(dest = path_cur$angle_record, orig = ess_angles)
}