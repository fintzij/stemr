#' Initialize the ODE path
#'
#' @param data matrix containing the dataset
#' @param parmat matrix with parameters, contants, time-varying pars and covars
#' @param param_blocks list of parameter blocks
#' @param tparam list of time-varying parameters
#' @param censusmat template matrix for the LNA path and incidence at the
#'   observation times
#' @param emitmat matrix in which to store the log-emission probabilities
#' @param stoich_matrix ODE stoichiometry matrix
#' @param proc_pointer external LNA pointer
#' @param set_pars_pointer pointer for setting the LNA parameters
#' @param times times at which the LNA should be evaluated
#' @param param_inds C++ column indices for parameters
#' @param const_inds C++ column indices for constants
#' @param tcovar_inds C++ column indices for time varying covariates
#' @param initdist_inds C++ column indices in the LNA parameter matrix for
#'   the initial state
#' @param param_update_inds logical vector indicating when to update the
#'   parameters
#' @param census_indices C++ row indices of LNA times when the path is to be
#'   censused
#' @param event_inds vector of column indices in the LNA path for which
#'   incidence will be computed.
#' @param measproc_indmat logical matrix for evaluating the measuement process
#' @param d_meas_pointer external pointer for the measurement process function
#' @param do_prevalence should prevalence be computed?
#' @param forcing_inds logical vector of indicating at which times in the
#'   time-varying covariance matrix a forcing is applied.
#' @param initialization_attempts number of initialization attempts
#' @param step_size initial step size for the ODE solver (adapted internally,
#'   but too large of an initial step can lead to failure in stiff systems).
#' @param par_init_fcn function for initializing the parameter values
#' @param param_vec vector for storing lna parameters when evaluating the
#'   measurement process
#'
#' @return ODE path
#' @export
initialize_ode <-
        function(dat,
                 parmat,
                 param_blocks,
                 tparam,
                 censusmat,
                 emitmat,
                 stoich_matrix,
                 proc_pointer,
                 set_pars_pointer,
                 times,
                 param_vec,
                 param_inds,
                 const_inds,
                 tcovar_inds,
                 initdist_inds,
                 param_update_inds,
                 census_indices,
                 event_inds,
                 measproc_indmat,
                 d_meas_pointer,
                 do_prevalence,
                 forcing_inds,
                 forcing_tcov_inds,
                 forcings_out,
                 forcing_transfers,
                 initialization_attempts,
                 step_size,
                 fixed_inits,
                 init_volumes_cur,
                 initdist_objects) {
              
              # get the initial state parameters
              init_state   <- parmat[1, initdist_inds + 1]
              data_log_lik <- NaN
              attempt      <- 0
              keep_going   <- TRUE
              flow_matrix  <- t(stoich_matrix)
              
              while(keep_going && (attempt <= initialization_attempts)) {
                    
                    try({
                          # integrate the ODEs
                          path_init <- 
                              integrate_odes(
                                  ode_times         = times,
                                  ode_pars          = parmat,
                                  ode_param_inds    = param_inds,
                                  ode_tcovar_inds   = tcovar_inds,
                                  init_start        = initdist_inds[1],
                                  param_update_inds = param_update_inds,
                                  stoich_matrix     = stoich_matrix,
                                  forcing_inds      = forcing_inds,
                                  forcing_tcov_inds = forcing_tcov_inds,
                                  forcings_out      = forcings_out,
                                  forcing_transfers = forcing_transfers,
                                  step_size         = step_size,
                                  ode_pointer       = proc_pointer,
                                  set_pars_pointer  = set_pars_pointer)
                          
                          path <- list(latent_path = path_init$incid_path)
                          
                          census_latent_path(
                              path                = path$latent_path,
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
                          
                          # compute the dat log likelihood
                          data_log_lik <- sum(emitmat[,-1][measproc_indmat])
                          if(is.nan(data_log_lik)) data_log_lik <- -Inf
                    }, silent = TRUE)
                    
                    # propose new parameter values and/or initial volumes
                    keep_going <- is.nan(data_log_lik) || data_log_lik == -Inf
                    attempt    <- attempt + 1
                    
                    # try new parameters
                    if(keep_going) {
                        
                        for(s in seq_along(initdist_objects)) {
                            
                            if(!initdist_objects[[s]]$fixed) {
                                
                                # N(0,1) draws
                                draw_normals(initdist_objects[[s]]$draws_cur)
                                
                                # map draws
                                copy_vec(dest = initdist_objects[[s]]$init_volumes,
                                         orig = c(initdist_objects[[s]]$comp_mean +
                                                      c(initdist_objects[[s]]$comp_sqrt_cov %*%
                                                            initdist_objects[[s]]$draws_cur)))
                                
                                while(any(initdist_objects[[s]]$init_volumes < 0) | 
                                      any(init_volumes_cur[initdist_objects[[s]]$comp_inds_R] > initdist_objects[[s]]$comp_size)) {
                                    
                                    # N(0,1) draws
                                    draw_normals(initdist_objects[[s]]$draws_cur)
                                    
                                    # map draws
                                    copy_vec(dest = initdist_objects[[s]]$init_volumes,
                                             orig = c(initdist_objects[[s]]$comp_mean +
                                                          c(initdist_objects[[s]]$comp_sqrt_cov %*%
                                                                initdist_objects[[s]]$draws_cur))) 
                                }
                            }
                        }
                        
                        # copy to the ode parameter matrix
                        insert_initdist(parmat = parmat,
                                        initdist_objects = initdist_objects, 
                                        prop = FALSE)
                                      
                        # draw new parameter values if called for
                        for(s in seq_along(param_blocks)) {
                            if(!is.null(param_blocks[[s]]$initializer)) {
                                
                                # initialize parameters
                                param_blocks[[s]]$pars_nat = 
                                    param_blocks[[s]]$initializer()
                                
                                param_blocks[[s]]$pars_est = 
                                    param_blocks[[s]]$priors$to_estimation_scale(
                                        param_blocks[[s]]$pars_nat)
                                
                                # calculate log prior
                                param_blocks[[s]]$log_pd = 
                                    param_blocks[[s]]$priors$logprior(
                                        param_blocks[[s]]$pars_est)
                            }
                        }
                        
                        # insert parameters into the parameter matrix
                        insert_params(parmat = parmat,
                                      param_blocks = param_blocks)
                                
                        if(!is.null(tparam)) {
                            for(s in seq_along(tparam)) {
                                
                                # sample new draws
                                draw_normals(tparam[[s]]$draws_cur)
                                tparam[[s]]$log_lik <- 
                                    sum(dnorm(tparam[[s]]$draws_cur, log = TRUE))
                                
                                # get values
                                insert_tparam(tcovar    = parmat,
                                              values    = 
                                                  tparam[[s]]$draws2par(
                                                      parameters = parmat[1,], 
                                                      draws = tparam[[s]]$draws_cur),
                                              col_ind   = tparam[[s]]$col_ind,
                                              tpar_inds = tparam[[s]]$tpar_inds)
                            }
                        }   
                    }
              }
              
              if(keep_going) {
                    
                    stop("Initialization failed. Try different initial parameter values.")
                    
              } else {
                    
                  path$data_log_lik <- data_log_lik # sum of log emission probabilities
                  return(list(path = path,
                              param_blocks = param_blocks,
                              initdist_objects = initdist_objects,
                              tparam = tparam))
              }
        }