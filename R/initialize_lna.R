#' Initialize the LNA path
#'
#' @param data matrix containing the dataset
#' @param lna_parameters parameters, contants, time-varying covariates at LNA
#'   times
#' @param tparam list of time-varying parameters
#' @param censusmat template matrix for the LNA path and incidence at the
#'   observation times
#' @param emitmat matrix in which to store the log-emission probabilities
#' @param stoich_matrix LNA stoichiometry matrix
#' @param lna_pointer external LNA pointer
#' @param lna_set_pars_pointer pointer for setting the LNA parameters
#' @param lna_times times at whicht eh LNA should be evaluated
#' @param lna_param_inds C++ column indices for parameters
#' @param lna_const_inds C++ column indices for constants
#' @param lna_tcovar_inds C++ column indices for time varying covariates
#' @param lna_initdist_inds C++ column indices in the LNA parameter matrix for
#'   the initial state
#' @param param_update_inds logical vector indicating when to update the
#'   parameters
#' @param census_indices C++ row indices of LNA times when the path is to be
#'   censused
#' @param lna_event_inds vector of column indices in the LNA path for which
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
#' @param par_init_fcn function for initializing the parameter values
#' @param ess_warmup number of elliptical slice sampling updates where
#'   likelihood is over indicators for monotonicity and non-negativity of LNA
#'   increments
#' @param lna_param_vec vector for storing lna parameters when evaluating the
#'   measurement process
#'
#' @return LNA path along with its stochastic perturbations
#' @export
initialize_lna <-
        function(data,
                 lna_parameters,
                 tparam,
                 censusmat,
                 emitmat,
                 stoich_matrix,
                 lna_pointer,
                 lna_set_pars_pointer,
                 lna_times,
                 lna_param_vec,
                 lna_param_inds,
                 lna_const_inds,
                 lna_tcovar_inds,
                 lna_initdist_inds,
                 param_update_inds,
                 census_indices,
                 lna_event_inds,
                 measproc_indmat,
                 d_meas_pointer,
                 do_prevalence,
                 forcing_inds,
                 forcing_matrix,
                 initialization_attempts,
                 step_size,
                 par_init_fcn = NULL,
                 ess_warmup) {

                # get the initial state parameters
                init_state   <- lna_parameters[1, lna_initdist_inds + 1]
                data_log_lik <- NaN
                attempt      <- 0
                keep_going   <- TRUE

                while(keep_going && (attempt <= initialization_attempts)) {
                        try({
                                # propose another LNA path
                                path_init <- propose_lna_approx(
                                        lna_times         = lna_times,
                                        lna_pars          = lna_parameters,
                                        init_start        = lna_initdist_inds[1],
                                        param_update_inds = param_update_inds,
                                        stoich_matrix     = stoich_matrix,
                                        forcing_inds      = forcing_inds,
                                        forcing_matrix    = forcing_matrix,
                                        max_attempts      = initialization_attempts,
                                        step_size         = step_size, 
                                        nsim              = 1, 
                                        ess_updates       = 1, 
                                        ess_warmup        = ess_warmup,
                                        lna_pointer       = lna_pointer,
                                        set_pars_pointer  = lna_set_pars_pointer
                                )
                                
                                path <- list(draws    = t(path_init$draws[,,1]),
                                             lna_path = path_init$incid_paths[,,1])

                                census_lna(
                                        path                = path$lna_path,
                                        census_path         = censusmat,
                                        census_inds         = census_indices,
                                        lna_event_inds      = lna_event_inds,
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
                                data_log_lik <- sum(emitmat[,-1][measproc_indmat])
                                if(is.nan(data_log_lik)) data_log_lik <- -Inf
                        }, silent = TRUE)

                        keep_going <- is.nan(data_log_lik) || data_log_lik == -Inf
                        attempt    <- attempt + 1

                        if(keep_going && !is.null(par_init_fcn)) {
                                pars2lnapars(lna_parameters, par_init_fcn())
                        }
                        
                        if(keep_going && !is.null(tparam)) {
                              
                              for(s in seq_along(tparam)) {
                                    
                                    # sample new draws
                                    draw_normals(tparam[[s]]$draws_cur)
                                    
                                    # get values
                                    insert_tparam(tcovar    = lna_parameters,
                                                  values    = tparam[[s]]$draws2par(parameters = lna_parameters[1,], draws = tparam[[s]]$draws_cur),
                                                  col_ind   = tparam[[s]]$col_ind,
                                                  tpar_inds = tparam[[s]]$tpar_inds)
                              }
                        }
                }

                if(keep_going) {

                        stop("Initialization failed. Try different initial parameter values.")

                } else {

                        path$data_log_lik <- data_log_lik # sum of log emission probabilities
                        return(path)
                }
        }