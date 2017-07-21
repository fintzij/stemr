#' Sample a new LNA path via elliptical slice sampling.
#'
#' @param path_cur list with the current LNA path along with its ODE paths
#' @param n_ess_updates number of elliptical slice sampling updates
#' @inheritParams initialize_lna
#'
#' @return list with an updated LNA path along with its stochastic
#'   perturbations, the observed data log-likelihood, and the lna
#'   log-likelihood, and a record of the number of elliptical slice sampling
#'   proposals
#' @export
update_lna_path <-
        function(path_cur, data, lna_parameters, pathmat_prop, censusmat, draws_prop, emitmat,
                 flow_matrix, stoich_matrix, lna_pointer, lna_set_pars_pointer, lna_times,
                 lna_params_temp, lna_param_inds, lna_const_inds, lna_tcovar_inds,
                 lna_initdist_inds, param_update_inds, incidence_codes,
                 census_incidence_codes, census_indices, measproc_indmat, obstime_inds,
                 d_meas_pointer, do_prevalence, n_ess_updates) {

        # vector for storing the number of steps in each ESS update
        ess_record <- rep(1, n_ess_updates)

        # get the initial state parameters and census the LNA path
        init_state <- lna_parameters[1, lna_initdist_inds + 1]

        # perform the ESS updates, retaining the last state
        for(k in seq_len(n_ess_updates)) {

                # sample a new set of stochastic perturbations
                perturbations <- matrix(rnorm(path_cur$draws),
                                     nrow = nrow(path_cur$draws),
                                     ncol = ncol(path_cur$draws))

                # choose a likelihood threshold
                threshold <- path_cur$data_log_lik + log(runif(1))

                # initial proposal, which also defines a bracket
                theta <- runif(1, 0, 2*pi)
                lower <- theta - 2*pi; upper <- theta

                # construct the first proposal
                copy_mat(draws_prop, cos(theta)*path_cur$draws + sin(theta)*perturbations)

                # initialize the data log likelihood for the proposed path
                data_log_lik_prop <- NULL

                # map the perturbations to an LNA path
                try({
                        map_draws_2_lna(pathmat           = pathmat_prop,
                                        draws             = draws_prop,
                                        lna_times         = lna_times,
                                        lna_pars          = lna_parameters,
                                        init_start        = lna_initdist_inds[1],
                                        param_update_inds = param_update_inds,
                                        stoich_matrix     = stoich_matrix,
                                        lna_pointer       = lna_pointer,
                                        set_pars_pointer  = lna_set_pars_pointer
                        )

                        census_lna(
                                path                = pathmat_prop,
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
                                lna_params_temp   = lna_params_temp,
                                lna_param_inds    = lna_param_inds,
                                lna_const_inds    = lna_const_inds,
                                lna_tcovar_inds   = lna_tcovar_inds,
                                param_update_inds = param_update_inds,
                                census_indices    = census_indices,
                                d_meas_ptr        = d_meas_pointer
                        )

                        # compute the data log likelihood
                        data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                        if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
                }, silent = TRUE)

                if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf

                # accept or reject the proposal
                while((data_log_lik_prop < threshold) && !isTRUE(all.equal(lower, upper))) {

                        # increment the number of ESS proposals for the current iteration
                        ess_record[k] <- ess_record[k] + 1

                        # shrink the bracket
                        if(theta < 0) {
                                lower <- theta
                        } else {
                                upper <- theta
                        }

                        # sample a new point
                        theta <- runif(1, lower, upper)

                        # construct the next proposal
                        copy_mat(draws_prop, cos(theta)*path_cur$draws + sin(theta)*perturbations)

                        # initialize the data log likelihood for the proposed path
                        data_log_lik_prop <- NULL

                        # map the perturbations to an LNA path
                        try({
                                map_draws_2_lna(pathmat           = pathmat_prop,
                                                draws             = draws_prop,
                                                lna_times         = lna_times,
                                                lna_pars          = lna_parameters,
                                                init_start        = lna_initdist_inds[1],
                                                param_update_inds = param_update_inds,
                                                stoich_matrix     = stoich_matrix,
                                                lna_pointer       = lna_pointer,
                                                set_pars_pointer  = lna_set_pars_pointer
                                )

                                census_lna(
                                        path                = pathmat_prop,
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
                                        lna_params_temp   = lna_params_temp,
                                        lna_param_inds    = lna_param_inds,
                                        lna_const_inds    = lna_const_inds,
                                        lna_tcovar_inds   = lna_tcovar_inds,
                                        param_update_inds = param_update_inds,
                                        census_indices    = census_indices,
                                        d_meas_ptr        = d_meas_pointer
                                )

                                # compute the data log likelihood
                                data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                                if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
                        }, silent = TRUE)

                        if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
                }

                # if the bracket width is not equal to zero, update the draws, path, and data log likelihood
                if(!isTRUE(all.equal(lower, upper))) {
                        # transfer the new path and residual path into the* sin(theta) path_prop list
                        copy_mat(path_cur$draws, draws_prop)
                        copy_mat(path_cur$lna_path, pathmat_prop)
                        copy_vec(path_cur$data_log_lik, data_log_lik_prop)
                }
        }

        copy_vec(path_cur$ess_record, ess_record)
}