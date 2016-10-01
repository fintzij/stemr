#' Initialize the latent LNA path
#'
#' @inheritParams propose_lna_path
#' @param incidence_codes vector of incidence codes
#' @param data matrix containing the data
#' @param measproc_indmat matrix indicating which compartments are measured
#' @param obstime_inds list of indices indicating at which rows each measured
#'   variable is recorded
#' @param d_meas_pointer external pointer to density function
#' @param parameters vector of model parameters
#' @param constants vector of constants
#' @param tcovar_censmat time-varying covariates at observation times
#' @param initialization_attempts maximum number of initialization attempts
#'
#' @return list containing an initial latent path along with the ODEs
#' @export
initialize_lna <- function(lna_parameters, stoich_matrix, lna_pointer, lna_times, initdist_parameters, fixed_inits, param_update_inds, drift_inds, resid_inds, diff_inds, log_scale, incidence_codes, data, measproc_indmat, obstime_inds, d_meas_pointer, parameters, constants, tcovar_censmat, initialization_attempts) {

        path <- propose_lna_path(parameters         = lna_parameters,
                                 stoich_matrix       = stoich_matrix,
                                 lna_pointer         = lna_pointer,
                                 times               = lna_times,
                                 initdist_parameters = initdist_parameters,
                                 fixed_inits         = fixed_inits,
                                 param_update_inds   = param_update_inds,
                                 drift_inds          = drift_inds,
                                 resid_inds          = resid_inds,
                                 diff_inds           = diff_inds,
                                 log_scale           = log_scale,
                                 incidence_codes     = incidence_codes)

        # matrix in which to store the emission probabilities
        emitmat <- cbind(data[,1,drop = F], matrix(0.0, nrow = nrow(measproc_indmat),
                                                   ncol = ncol(measproc_indmat), dimnames = list(NULL, colnames(measproc_indmat))))

        # does the path need to be censused at observation times (e.g. if there
        # are time-varying covariates that change at times other than
        # observation times)
        do_census    <- !identical(lna_times, data[,1])
        do_incidence <- !is.null(incidence_codes)

        if(do_census) {
                statemat <- path$path[match(data[,1], path$path[,1]),]
        } else {
                statemat <- path$path
        }

        if(log_scale) {
                statemat[,-1] <- exp(statemat[,-1])
        }

        if(do_incidence) {
                compute_incidence(statemat, incidence_codes+1, row_inds = obstime_inds)
        }

        # evaluate the density
        evaluate_d_measure(emitmat          = emitmat,
                           obsmat           = data,
                           statemat         = statemat,
                           measproc_indmat  = measproc_indmat,
                           parameters       = parameters,
                           constants        = constants,
                           tcovar_censusmat = tcovar_censmat,
                           d_meas_ptr       = d_meas_pointer)

        attempt    <- 1
        keep_going <- (any(is.nan(emitmat)) || any(emitmat == -Inf))

        while(attempt <= initialization_attempts && keep_going) {

                path <- propose_lna_path(parameters         = lna_parameters,
                                         stoich_matrix       = stoich_matrix,
                                         lna_pointer         = lna_pointer,
                                         times               = lna_times,
                                         initdist_parameters = initdist_parameters,
                                         fixed_inits         = fixed_inits,
                                         param_update_inds   = param_update_inds,
                                         drift_inds          = drift_inds,
                                         resid_inds          = resid_inds,
                                         diff_inds           = diff_inds,
                                         log_scale           = log_scale,
                                         incidence_codes     = incidence_codes)

                if(do_census) {
                        statemat <- path$path[match(data[,1], path$path[,1]),]
                } else {
                        statemat <- path$path
                }


                if(log_scale) {
                        statemat[,-1] <- exp(statemat[,-1])
                }

                if(do_incidence) {
                        compute_incidence(statemat, incidence_codes+1, row_inds = obstime_inds)
                }

                # evaluate the density
                evaluate_d_measure(emitmat          = emitmat,
                                   obsmat           = data,
                                   statemat         = statemat,
                                   measproc_indmat  = measproc_indmat,
                                   parameters       = parameters,
                                   constants        = constants,
                                   tcovar_censusmat = tcovar_censmat,
                                   d_meas_ptr       = d_meas_pointer)

                attempt <- attempt + 1
                keep_going <- any(is.nan(emitmat)) || any(emitmat == -Inf)
        }

        if(keep_going) {
                stop("Initialization failed. Try different initial parameter values.")
        } else {
                return(list(path = path, data_log_lik = sum(emitmat[,-1][measproc_indmat])))
        }
}