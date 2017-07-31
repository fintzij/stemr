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
#' @param incidence_codes C++ indices for LNA count process compartments for
#'   which incidence is to be computed
#' @param census_incidence_codes C++ indices for columns in the census matrix
#'   containing incidence data
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
                 incidence_codes,
                 census_incidence_codes,
                 census_indices,
                 measproc_indmat,
                 obstime_inds,
                 d_meas_pointer,
                 do_prevalence,
                 initialization_attempts) {

        # get the initial state parameters
        init_state <- lna_parameters[1, lna_initdist_inds + 1]

        data_log_lik <- NULL

        try({
                # propose another LNA path
                path <- propose_lna(
                        lna_times            = lna_times,
                        lna_pars             = lna_parameters,
                        init_start           = stem_object$dynamics$lna_initdist_inds[1],
                        param_update_inds    = param_update_inds,
                        stoich_matrix        = stoich_matrix,
                        lna_pointer          = lna_pointer,
                        set_pars_pointer     = lna_set_pars_pointer
                )

                census_lna(
                        path                = path$lna_path,
                        census_path         = censusmat,
                        census_inds         = census_indices,
                        flow_matrix_lna     = t(stoich_matrix),
                        do_prevalence       = do_prevalence,
                        init_state          = init_state,
                        incidence_codes_lna = incidence_codes
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

        attempt <- 1
        keep_going <- is.nan(data_log_lik) || data_log_lik == -Inf

        while(attempt <= initialization_attempts && keep_going) {

                data_log_lik <- NULL

                try({
                        # propose another LNA path
                        path <- propose_lna(
                                lna_times            = lna_times,
                                lna_pars             = lna_parameters,
                                init_start           = stem_object$dynamics$lna_initdist_inds[1],
                                param_update_inds    = param_update_inds,
                                stoich_matrix        = stoich_matrix,
                                lna_pointer          = lna_pointer,
                                set_pars_pointer     = lna_set_pars_pointer
                        )

                        census_lna(
                                path                = path$lna_path,
                                census_path         = censusmat,
                                census_inds         = census_indices,
                                flow_matrix_lna     = t(stoich_matrix),
                                do_prevalence       = do_prevalence,
                                init_state          = init_state,
                                incidence_codes_lna = incidence_codes
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

                attempt <- attempt + 1
                keep_going <- is.nan(data_log_lik) || data_log_lik == -Inf
        }

        if(keep_going) {
                stop("Initialization failed. Try different initial parameter values.")
        } else {
                path$data_log_lik <- data_log_lik # sum of log emission probabilities
                return(path)
        }
}