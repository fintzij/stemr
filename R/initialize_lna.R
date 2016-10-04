#' Initialize the latent LNA path
#'
#' @inheritParams propose_lna
#' @param incidence_codes_lna vector of incidence codes
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
initialize_lna <- function(lna_parameters, flow_matrix, lna_pointer, lna_set_pars_pointer, lna_times, fixed_inits, param_update_inds, incidence_codes, data, measproc_indmat, obstime_inds, d_meas_pointer, parameters, constants, tcovar_censmat, initialization_attempts) {

        path <- propose_lna(lna_times         = lna_times,
                            lna_pars          = lna_parameters,
                            param_update_inds = param_update_inds,
                            flow_matrix       = flow_matrix,
                            lna_pointer       = lna_pointer,
                            set_pars_pointer  = lna_set_pars_pointer)

        # matrix in which to store the emission probabilities
        emitmat <- cbind(data[,1,drop = F],
                         matrix(0.0, nrow = nrow(measproc_indmat), ncol = ncol(measproc_indmat),
                                dimnames = list(NULL, colnames(measproc_indmat))))

        # does the path need to be censused at observation times (e.g. if there
        # are time-varying covariates that change at times other than
        # observation times)
        do_census    <- !identical(lna_times, data[,1])
        do_incidence <- !is.null(incidence_codes)

        # get the initial state parameters
        init_state   <- lna_parameters[1, grepl("_0", colnames(lna_parameters))]

        # generate the statemat
        if(do_census) {
                statemat <- path$lna_path[match(data[,1], path$lna_path[,1]),]
        } else {
                statemat <- path$lna_path
        }

        # either compute the incidence, or compute the compartment counts if the data are prevalence counts
        if(do_incidence) {
                # compute the incidence
                compute_incidence(statemat, incidence_codes+1, row_inds = obstime_inds)

                # evaluate the density of the incidence counts
                evaluate_d_measure(emitmat          = emitmat,
                                   obsmat           = data,
                                   statemat         = statemat,
                                   measproc_indmat  = measproc_indmat,
                                   parameters       = parameters,
                                   constants        = constants,
                                   tcovar_censusmat = tcovar_censmat,
                                   d_meas_ptr       = d_meas_pointer)

        } else {
                # convert the event counts to the natrual state space
                pathmat <- convert_lna(statemat, flow_matrix, init_state)

                # evaluate the density of the prevalence counts
                evaluate_d_measure(emitmat          = emitmat,
                                   obsmat           = data,
                                   statemat         = pathmat,
                                   measproc_indmat  = measproc_indmat,
                                   parameters       = parameters,
                                   constants        = constants,
                                   tcovar_censusmat = tcovar_censmat,
                                   d_meas_ptr       = d_meas_pointer)

        }

        attempt    <- 1
        keep_going <- (any(is.nan(emitmat)) || any(emitmat == -Inf))

        while(attempt <= initialization_attempts && keep_going) {

                path <- propose_lna(lna_times         = lna_times,
                                    lna_pars          = lna_parameters,
                                    param_update_inds = param_update_inds,
                                    flow_matrix       = flow_matrix,
                                    lna_pointer       = lna_pointer,
                                    set_pars_pointer  = lna_set_pars_pointer)

                if(do_census) {
                        statemat <- path$lna_path[match(data[,1], path$lna_path[,1]),]
                } else {
                        statemat <- path$lna_path
                }

                # either compute the incidence, or compute the compartment counts if the data are prevalence counts
                if(do_incidence) {
                        # compute the incidence
                        compute_incidence(statemat, incidence_codes+1, row_inds = obstime_inds)

                        # evaluate the density of the incidence counts
                        evaluate_d_measure(emitmat          = emitmat,
                                           obsmat           = data,
                                           statemat         = statemat,
                                           measproc_indmat  = measproc_indmat,
                                           parameters       = parameters,
                                           constants        = constants,
                                           tcovar_censusmat = tcovar_censmat,
                                           d_meas_ptr       = d_meas_pointer)

                } else {
                        # convert the event counts to the natrual state space
                        convert_lna2(statemat, flow_matrix, init_state, pathmat)

                        # evaluate the density of the prevalence counts
                        evaluate_d_measure(emitmat          = emitmat,
                                           obsmat           = data,
                                           statemat         = pathmat,
                                           measproc_indmat  = measproc_indmat,
                                           parameters       = parameters,
                                           constants        = constants,
                                           tcovar_censusmat = tcovar_censmat,
                                           d_meas_ptr       = d_meas_pointer)

                }

                attempt <- attempt + 1
                keep_going <- any(is.nan(emitmat)) || any(emitmat == -Inf)
        }

        if(keep_going) {
                stop("Initialization failed. Try different initial parameter values.")
        } else {
                path$data_log_lik = sum(emitmat[,-1][measproc_indmat])
                return(path)
        }
}