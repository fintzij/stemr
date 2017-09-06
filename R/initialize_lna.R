#' Initialize the LNA path
#'
#' @param data matrix containing the dataset
#' @param lna_parameters parameters, contants, time-varying covariates at LNA
#'   times
#' @param censusmat template matrix for the LNA path and incidence at the
#'   observation times
#' @param emitmat matrix in which to store the log-emission probabilities
#' @param stoich_matrix LNA stoichiometry matrix
#' @param lna_pointer external LNA pointer
#' @param lna_set_pars_pointer pointer for setting the LNA parameters
#' @param lna_times times at whicht eh LNA should be evaluated
#' @param lna_initdist_inds C++ column indices in the LNA parameter matrix for
#'   the initial state
#' @param param_update_inds logical vector indicating when to update the
#'   parameters
#' @param census_indices C++ row indices of LNA times when the path is to be
#'   censused
#' @param measproc_indmat logical matrix for evaluating the measuement process
#' @param obstime_inds list of C++ indices indicating when each measurement
#'   process is evaluated
#' @param d_meas_pointer external pointer for the measurement process function
#' @param parameters vector of model parameters
#' @param constants vector of model constants
#' @param tcovar_censmat time-varying covariates at LNA times
#' @param do_prevalence should prevalence be computed?
#' @param forcing_inds logical vector of indicating at which times in the
#'   time-varying covariance matrix a forcing is applied.
#' @param forcing_matrix matrix containing the forcings.
#' @param initialization_attempts number of initialization attempts
#'
#' @return LNA path along with its stochastic perturbations
#' @export
initialize_lna <-
        function(data,
                 lna_parameters,
                 censusmat,
                 emitmat,
                 stoich_matrix,
                 lna_pointer,
                 lna_set_pars_pointer,
                 lna_times,
                 lna_param_inds,
                 lna_const_inds,
                 lna_tcovar_inds,
                 lna_initdist_inds,
                 param_update_inds,
                 census_indices,
                 measproc_indmat,
                 obstime_inds,
                 d_meas_pointer,
                 do_prevalence,
                 forcing_inds,
                 forcing_matrix,
                 initialization_attempts) {

        # get the initial state parameters
        init_state <- lna_parameters[1, lna_initdist_inds + 1]

        data_log_lik <- NaN
        attempt      <- 0
        keep_going   <- TRUE

        while(keep_going && (attempt <= initialization_attempts)) {
                try({
                        # propose another LNA path
                        path <- propose_lna(
                                lna_times         = lna_times,
                                lna_pars          = lna_parameters,
                                init_start        = stem_object$dynamics$lna_initdist_inds[1],
                                param_update_inds = param_update_inds,
                                stoich_matrix     = stoich_matrix,
                                forcing_inds      = forcing_inds,
                                forcing_matrix    = forcing_matrix,
                                step_size         = stem_object$dynamics$dynamics_args$step_size,
                                lna_pointer       = stem_object$dynamics$lna_pointers$lna_ptr,
                                set_pars_pointer  = stem_object$dynamics$lna_pointers$set_lna_params_ptr
                        )

                        path$prev_path <- NULL

                        census_lna(
                                path                = path$lna_path,
                                census_path         = censusmat,
                                census_inds         = census_indices,
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
                                d_meas_ptr        = d_meas_pointer
                        )

                        # compute the data log likelihood
                        data_log_lik <- sum(emitmat[,-1][measproc_indmat])
                        if(is.nan(data_log_lik)) data_log_lik <- -Inf
                }, silent = TRUE)

                keep_going <- is.nan(data_log_lik) || data_log_lik == -Inf
                attempt    <- attempt + 1
        }

        if(keep_going) {

                stop("Initialization failed. Try different initial parameter values.")

        } else {

                path$data_log_lik <- data_log_lik # sum of log emission probabilities
                return(path)
        }
}