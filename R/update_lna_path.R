#' Sample a new LNA path via elliptical slice sampling.
#'
#' @param path_cur list with the current LNA path along with its moments
#' @inheritParams initialize_lna
#' @param do_census Does the path need to be censused?
#' @param do_incidence Does incidence need to be computed?
#'
#' @return list with an updated LNA path along with its ODEs
#' @export
update_lna_path <- function(path_cur, lna_parameters, stoich_matrix, lna_ess_pointer, lna_times, initdist_parameters, fixed_inits, param_update_inds, drift_inds, resid_inds, log_scale, incidence_codes, data, measproc_indmat, obstime_inds, d_meas_pointer, parameters, constants, tcovar_censmat, do_census, do_incidence) {

        # matrix in which to store the emission probabilities
        emitmat     <- matrix(0.0, nrow = nrow(measproc_indmat), ncol = ncol(measproc_indmat)+1,
                              dimnames = list(NULL, c("time", colnames(measproc_indmat))))
        emitmat[,1] <- data[,1]
        pathmat     <- matrix(0.0, nrow = nrow(path_cur$path$path), ncol = ncol(path_cur$path$path))
        pathmat[,1] <- lna_times

        # propose a new path
        path_prop <- propose_lna_path_ess(path_cur            = path_cur,
                                          parameters          = lna_parameters,
                                          stoich_matrix       = stoich_matrix,
                                          lna_ess_pointer     = lna_ess_pointer,
                                          times               = lna_times,
                                          initdist_parameters = initdist_parameters,
                                          fixed_inits         = fixed_inits,
                                          param_update_inds   = param_update_inds,
                                          drift_inds          = drift_inds,
                                          resid_inds          = resid_inds,
                                          log_scale           = log_scale,
                                          incidence_codes     = incidence_codes)

        # choose a likelihood threshold
        threshold <- path_cur$data_log_lik + log(runif(1))

        # initial proposal, which also defines a bracket
        theta <- runif(1, 0, 2*pi)
        lower <- theta - 2*pi; upper <- theta

        # construct the first proposal
        pathmat[, -1] <- path_prop$drift_process + cos(theta)*path_cur$path$residual_path + sin(theta)*path_prop$residual_path

        # does the path_new need to be censused at observation times (e.g. if there
        # are time-varying covariates that change at times other than
        # observation times)

        if(do_census) {
                statemat <- pathmat[match(data[,1], pathmat[,1]),]
        } else {
                statemat <- pathmat
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

        # compute the data log likelihood
        data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
        if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf

        # accept or reject the proposal
        while(data_log_lik_prop < threshold) {

                # shrink the bracket
                if(theta < 0) {
                        lower <- theta
                } else {
                        upper <- theta
                }

                # sample a new point
                theta <- runif(1, lower, upper)

                # construct the next proposal
                pathmat[, -1] <- path_prop$drift_process + cos(theta)*path_cur$path$residual_path + sin(theta)*path_prop$residual_path

                # does the path_new need to be censused at observation times (e.g. if there
                # are time-varying covariates that change at times other than
                # observation times)

                if(do_census) {
                        statemat <- pathmat[match(data[,1], pathmat[,1]),]
                } else {
                        statemat <- pathmat
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

                # compute the data log likelihood
                data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
        }

        return(list(path = list(path = pathmat,
                                residual_path = resmat,
                                drift_process = path_cur$path$drift_process,
                                residual_process = resmat,
                                diffusion_process=path_cur$path$diffusion_process),
                    data_log_lik = data_log_lik_prop))
}