#' Initialize the LNA path
#'
#' @param data matrix containing the dataset
#' @param ode_parameters parameters, contants, time-varying covariates at ode
#'   times
#' @param censusmat template matrix for the ode path and incidence at the
#'   observation times
#' @param emitmat matrix in which to store the log-emission probabilities
#' @param stoich_matrix ode stoichiometry matrix
#' @param ode_pointer external ode pointer
#' @param ode_set_pars_pointer pointer for setting the ode parameters
#' @param ode_times times at whicht eh ode should be evaluated
#' @param ode_param_inds C++ column indices for parameters
#' @param ode_const_inds C++ column indices for constants
#' @param ode_tcovar_inds C++ column indices for time varying covariates
#' @param ode_initdist_inds C++ column indices in the ode parameter matrix for
#'   the initial state
#' @param param_update_inds logical vector indicating when to update the
#'   parameters
#' @param census_indices C++ row indices of ode times when the path is to be
#'   censused
#' @param ode_event_inds vector of column indices in the ode path for which
#'   incidence will be computed.
#' @param measproc_indmat logical matrix for evaluating the measuement process
#' @param d_meas_pointer external pointer for the measurement process function
#' @param do_prevalence should prevalence be computed?
#' @param forcing_inds logical vector of indicating at which times in the
#'   time-varying covariance matrix a forcing is applied.
#' @param forcing_matrix matrix containing the forcings.
#' @param initialization_attempts number of initialization attempts
#' @param step_size initial step size for the ODE solver (adapted internally,
#'   but too large of an initial step can lead to failure in stiff systems).
#'
#' @return ODE path
#' @export
initialize_ode <-
        function(data,
                 ode_parameters,
                 tparam,
                 censusmat,
                 emitmat,
                 stoich_matrix,
                 ode_pointer,
                 ode_set_pars_pointer,
                 ode_times,
                 ode_param_vec,
                 ode_param_inds,
                 ode_const_inds,
                 ode_tcovar_inds,
                 ode_initdist_inds,
                 param_update_inds,
                 census_indices,
                 ode_event_inds,
                 measproc_indmat,
                 d_meas_pointer,
                 do_prevalence,
                 forcing_inds,
                 forcing_matrix,
                 initialization_attempts,
                 step_size,
                 par_init_fcn = NULL) {

        # get the initial state parameters
        init_state   <- ode_parameters[1, ode_initdist_inds + 1]
        data_log_lik <- NaN
        attempt      <- 0
        keep_going   <- TRUE

        while(keep_going && (attempt <= initialization_attempts)) {
              try({
                    # integrate the ODEs
                    path <- integrate_odes(ode_times         = ode_times,
                                           ode_pars          = ode_parameters,
                                           init_start        = ode_initdist_inds[1],
                                           param_update_inds = param_update_inds,
                                           stoich_matrix     = stoich_matrix,
                                           forcing_inds      = forcing_inds,
                                           forcing_matrix    = forcing_matrix,
                                           step_size         = step_size,
                                           ode_pointer       = ode_pointer,
                                           set_pars_pointer  = ode_set_pars_pointer
                    )
                    
                    path$prev_path <- NULL
                    names(path) <- "ode_path"
                    
                    census_lna(
                          path                = path$ode_path,
                          census_path         = censusmat,
                          census_inds         = census_indices,
                          lna_event_inds      = ode_event_inds,
                          flow_matrix_lna     = t(stoich_matrix),
                          do_prevalence       = do_prevalence,
                          init_state          = init_state,
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
                    data_log_lik <- sum(emitmat[,-1][measproc_indmat])
                    if(is.nan(data_log_lik)) data_log_lik <- -Inf
              }, silent = TRUE)
              
              if(is.null(par_init_fcn)) {
                    keep_going <- FALSE
                    
              } else {
                    keep_going <- is.nan(data_log_lik) || data_log_lik == -Inf
                    attempt    <- attempt + 1
                    
                    # try new parameters
                    pars2lnapars(ode_parameters, par_init_fcn())
                    if(!is.null(tparam)) {
                          
                          for(s in seq_along(tparam)) {
                                
                                # sample new draws
                                draw_normals(tparam[[s]]$draws_cur)
                                
                                # get values
                                insert_tparam(tcovar    = ode_parameters,
                                              values    = tparam[[s]]$draws2par(parameters = ode_parameters[1,], draws = tparam[[s]]$draws_cur),
                                              col_ind   = tparam[[s]]$col_ind,
                                              tpar_inds = tparam[[s]]$tpar_inds)
                          }
                    }
              }
        }
        
        if(data_log_lik == -Inf) {

                stop("Initialization failed. Try different initial parameter values.")

        } else {

                path$data_log_lik <- data_log_lik # sum of log emission probabilities
                return(path)
        }
}