#' Sample a new LNA path via elliptical slice sampling.
#'
#' @param path_cur list with the current LNA path along with its ODE paths
#' @param pathmat matrix in which to store the proposed LNA paths
#' @inheritParams initialize_lna
#'
#' @return list with an updated LNA path along with its ODEs, the observed data
#'   log-likelihood, and the lna log-likelihood.
#' @export
update_lna_path <- function(path_cur, data, lna_parameters, pathmat, censusmat, emitmat, flow_matrix, lna_pointer_ess, lna_ess_set_pars_ptr, lna_times, lna_initdist_inds, param_update_inds, incidence_codes, census_incidence_codes, census_indices, measproc_indmat, obstime_inds, d_meas_pointer, parameters, constants, tcovar_censmat, do_prevalence, do_incidence) {

        # propose a new path
        path_prop <- propose_lna_ess(
                path_cur          = path_cur,
                lna_times         = lna_times,
                lna_pars          = lna_parameters,
                param_update_inds = param_update_inds,
                flow_matrix       = flow_matrix,
                lna_pointer_ess   = lna_pointer_ess,
                lna_ess_set_pars_ptr = lna_ess_set_pars_ptr
        )

        # choose a likelihood threshold
        threshold <- path_cur$data_log_lik + log(runif(1))

        # initial proposal, which also defines a bracket
        theta <- runif(1, 0, 2*pi)
        lower <- theta - 2*pi; upper <- theta

        # get the initial state parameters
        init_state   <- lna_parameters[1, lna_initdist_inds + 1]

        # construct the first proposal
        pathmat[, -1] <- path_cur$drift + cos(theta)*path_cur$res_path[,-1] + sin(theta)*path_prop$res_path[,-1]

        # census the LNA
        census_lna(
                path                = pathmat,
                census_path         = censusmat,
                census_inds         = census_indices,
                flow_matrix_lna     = flow_matrix,
                do_prevalence       = do_prevalence,
                init_state          = init_state,
                incidence_codes_lna = incidence_codes
        )

        # either compute the incidence, or compute the compartment counts if the data are prevalence counts
        if(do_incidence) {
                # compute the incidence
                compute_incidence(censusmat = censusmat,
                                  col_inds  = census_incidence_codes,
                                  row_inds  = obstime_inds)
        }

        # evaluate the density of the incidence counts
        evaluate_d_measure(
                emitmat          = emitmat,
                obsmat           = data,
                statemat         = censusmat,
                measproc_indmat  = measproc_indmat,
                parameters       = parameters,
                constants        = constants,
                tcovar_censusmat = tcovar_censmat,
                d_meas_ptr       = d_meas_pointer
        )

        # compute the data log likelihood
        data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
        if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf

        # accept or reject the proposal
        while((data_log_lik_prop < threshold) && !isTRUE(all.equal(lower, upper))) {

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

                # census the LNA
                census_lna(
                        path                = pathmat,
                        census_path         = censusmat,
                        census_inds         = census_indices,
                        flow_matrix_lna     = flow_matrix,
                        do_prevalence       = do_prevalence,
                        init_state          = init_state,
                        incidence_codes_lna = incidence_codes
                )

                # either compute the incidence, or compute the compartment counts if the data are prevalence counts
                if(do_incidence) {
                        # compute the incidence
                        compute_incidence(censusmat = censusmat,
                                          col_inds  = census_incidence_codes,
                                          row_inds  = obstime_inds)
                }

                # evaluate the density of the incidence counts
                evaluate_d_measure(
                        emitmat          = emitmat,
                        obsmat           = data,
                        statemat         = censusmat,
                        measproc_indmat  = measproc_indmat,
                        parameters       = parameters,
                        constants        = constants,
                        tcovar_censusmat = tcovar_censmat,
                        d_meas_ptr       = d_meas_pointer

                )

                # compute the data log likelihood
                data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
        }

        # transfer the new path and residual path into the path_prop list
        path_prop$lna_path      <- pathmat
        path_prop$res_path[,-1] <- pathmat[,-1] - path_cur$drift
        path_prop$data_log_lik  <- data_log_lik_prop

        # compute the new LNA density
        path_prop <-
                lna_density2(
                        path_prop,
                        lna_times,
                        lna_parameters,
                        param_update_inds,
                        flow_matrix,
                        lna_pointer_ess,
                        lna_ess_set_pars_ptr
                )

        return(path_prop)
}