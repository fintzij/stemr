#' Compute the LNA log likelihood and data log likelihood.
#'
#' @param path list with the current LNA path along with its ODE paths to be
#'   reintegrated
#' @param params_prop_nat vector of proposed parameters on their natural scale
#' @param lna_params_prop matrix with proposed LNA parameters
#' @inheritParams initialize_lna
#'
#' @return list with an updated LNA path along with its ODEs, the observed data
#'   log-likelihood, and the lna log-likelihood.
#' @export
compute_lna_logliks <- function(path, data, params_prop_nat, lna_params_prop, censusmat, emitmat, flow_matrix, lna_times, lna_initdist_inds, param_update_inds, incidence_codes, census_incidence_codes, census_indices, measproc_indmat, obstime_inds, constants, tcovar_censmat, do_prevalence, do_incidence, lna_pointer, lna_set_pars_pointer, d_meas_pointer) {

        # Reintegrate the LNA ODEs and compute the LNA log likelihood
        path <- lna_density(
                path              = path,
                lna_times         = lna_times,
                lna_pars          = lna_params_prop,
                param_update_inds = param_update_inds,
                flow_matrix       = flow_matrix,
                lna_pointer       = lna_pointer,
                set_pars_pointer  = lna_set_pars_pointer)


        # get the initial state parameters
        init_state <- lna_params_prop[1, lna_initdist_inds + 1]

        # census the LNA
        census_lna(
                path                = path$lna_path,
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
                parameters       = params_prop_nat,
                constants        = constants,
                tcovar_censusmat = tcovar_censmat,
                d_meas_ptr       = d_meas_pointer
        )

        # compute the data log likelihood
        data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
        if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf

        # transfer the new path and residual path into the path list
        path$data_log_lik  <- data_log_lik_prop

        return(path)
}