#' Sample a new LNA path via elliptical slice sampling.
#'
#' @param path_cur list with the current LNA path along with its moments
#' @inheritParams initialize_lna
#' @param do_census Does the path need to be censused?
#' @param do_incidence Does incidence need to be computed?
#'
#' @return list with an updated LNA path along with its ODEs, the observed data
#'   log-likelihood, and the lna log-likelihood.
#' @export
update_lna_path <- function(path_cur, lna_parameters, flow_matrix, lna_ess_pointer, lna_times, lna_initdist_inds, param_update_inds, incidence_codes, data, measproc_indmat, obstime_inds, d_meas_pointer, parameters, constants, tcovar_censmat, do_census, do_incidence) {

        # matrix in which to store the emission probabilities
        emitmat     <- matrix(0.0, nrow = nrow(measproc_indmat), ncol = ncol(measproc_indmat)+1,
                              dimnames = list(NULL, c("time", colnames(measproc_indmat))))
        emitmat[,1] <- data[,1]
        pathmat     <- matrix(0.0, nrow = nrow(path_cur$lna_path), ncol = ncol(path_cur$lna_path))
        pathmat[,1] <- lna_times

        # propose a new path
        path_prop <- propose_lna_ess(path_cur          = path_cur,
                                     lna_times         = lna_times,
                                     lna_pars          = lna_parameters,
                                     param_update_inds = param_update_inds,
                                     flow_matrix       = flow_matrix,
                                     lna_pointer_ess   = lna_pointer_ess,
                                     set_pars_pointer  = lna_set_pars_pointer)

        # choose a likelihood threshold
        threshold <- path_cur$data_log_lik + log(runif(1))

        # initial proposal, which also defines a bracket
        theta <- runif(1, 0, 2*pi)
        lower <- theta - 2*pi; upper <- theta

        # get the initial state parameters
        init_state   <- lna_parameters[1, lna_initdist_inds + 1]

        # construct the first proposal
        pathmat[, -1] <- path_cur$drift + cos(theta)*path_cur$res_path[,-1] + sin(theta)*path_prop$res_path[,-1]

        # does the path_new need to be censused at observation times (e.g. if there
        # are time-varying covariates that change at times other than
        # observation times)
        if(do_census) {
                statemat <- pathmat[match(data[,1], pathmat[,1]),]
        } else {
                statemat <- pathmat
        }

        if(do_incidence) {
                # compute the incidence
                compute_incidence(statemat, incidence_codes+1, row_inds = obstime_inds)

                # evaluate the density
                evaluate_d_measure(emitmat          = emitmat,
                                   obsmat           = data,
                                   statemat         = statemat,
                                   measproc_indmat  = measproc_indmat,
                                   parameters       = parameters,
                                   constants        = constants,
                                   tcovar_censusmat = tcovar_censmat,
                                   d_meas_ptr       = d_meas_pointer)
        } else {
                # transform the path onto its natural state space
                natural_path <- convert_lna(statemat, flow_matrix, init_state)

                # evaluate the density
                evaluate_d_measure(emitmat          = emitmat,
                                   obsmat           = data,
                                   statemat         = natural_path,
                                   measproc_indmat  = measproc_indmat,
                                   parameters       = parameters,
                                   constants        = constants,
                                   tcovar_censusmat = tcovar_censmat,
                                   d_meas_ptr       = d_meas_pointer)
        }

        # compute the data log likelihood
        data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
        if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf

        # accept or reject the proposal
        while((data_log_lik_prop < threshold) && (abs(upper - lower)>.Machine$double.eps)) {

                # shrink the bracket
                if(theta < 0) {
                        lower <- theta
                } else {
                        upper <- theta
                }

                # sample a new point
                theta <- runif(1, lower, upper)

                # construct the next proposal
                pathmat[, -1] <- path_cur$drift + cos(theta)*path_cur$res_path[,-1] + sin(theta)*path_prop$res_path[,-1]

                # census the path if necessary
                if(do_census) {
                        statemat <- pathmat[match(data[,1], pathmat[,1]),]
                } else {
                        statemat <- pathmat
                }

                if(do_incidence) {
                        # compute the incidence
                        compute_incidence(statemat, incidence_codes+1, row_inds = obstime_inds)

                        # evaluate the density
                        evaluate_d_measure(emitmat          = emitmat,
                                           obsmat           = data,
                                           statemat         = statemat,
                                           measproc_indmat  = measproc_indmat,
                                           parameters       = parameters,
                                           constants        = constants,
                                           tcovar_censusmat = tcovar_censmat,
                                           d_meas_ptr       = d_meas_pointer)
                } else {
                        # transform the path onto its natural state space
                        natural_path <- convert_lna(statemat, flow_matrix, init_state)

                        # evaluate the density
                        evaluate_d_measure(emitmat          = emitmat,
                                           obsmat           = data,
                                           statemat         = natural_path,
                                           measproc_indmat  = measproc_indmat,
                                           parameters       = parameters,
                                           constants        = constants,
                                           tcovar_censusmat = tcovar_censmat,
                                           d_meas_ptr       = d_meas_pointer)
                }

                # compute the data log likelihood
                data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
        }

        # transfer the new path and residual path into the path_prop list
        path_prop$lna_path      <- pathmat
        path_prop$res_path[,-1] <- pathmat[,-1] - path_cur$drift
        path_prop$data_log_lik  <- data_log_lik_prop

        # compute the new LNA density
        path_prop <- lna_density2(path_prop, lna_times, lna_parameters,
                                  param_update_inds, flow_matrix, lna_pointer_ess, lna_set_pars_pointer)

        return(path_prop)
}