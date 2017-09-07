#' Sample a new LNA path via elliptical slice sampling.
#'
#' @param path_cur list with the current LNA path along with its ODE paths
#' @param n_ess_updates number of elliptical slice sampling updates
#' @param svd_sqrt,svd_d,svd_U,svd_V objects for computing the SVD of LNA
#'   diffusion matrics
#' @inheritParams initialize_lna
#'
#' @return list with an updated LNA path along with its stochastic
#'   perturbations, the observed data log-likelihood, and the lna
#'   log-likelihood, and a record of the number of elliptical slice sampling
#'   proposals
#' @export
update_lna_path <-
        function(path_cur,
                 data,
                 lna_parameters,
                 pathmat_prop,
                 censusmat,
                 draws_prop,
                 emitmat,
                 flow_matrix,
                 stoich_matrix,
                 lna_times,
                 forcing_inds,
                 forcing_matrix,
                 lna_param_inds,
                 lna_const_inds,
                 lna_tcovar_inds,
                 lna_initdist_inds,
                 param_update_inds,
                 census_indices,
                 lna_event_inds,
                 measproc_indmat,
                 svd_sqrt,
                 svd_d,
                 svd_U,
                 svd_V,
                 lna_pointer,
                 lna_set_pars_pointer,
                 d_meas_pointer,
                 do_prevalence,
                 n_ess_updates,
                 ess_schedule,
                 randomize_schedule,
                 step_size) {


        # vector for storing the number of steps in each ESS update
        ess_record <- matrix(1, nrow = n_ess_updates, ncol = length(ess_schedule[[1]]))

        # get the initial state parameters and census the LNA path
        init_state <- lna_parameters[1, lna_initdist_inds + 1]

        # perform the ESS updates, retaining the last state
        for(k in seq_len(n_ess_updates)) {

                if(randomize_schedule) {
                        ess_order <- sample.int(length(ess_schedule[[1]]), replace = FALSE)
                } else {
                        ess_order <- seq_along(ess_schedule[[1]])
                }

                for(j in ess_order) {

                        # sample a new set of stochastic perturbations
                        perturbations <- matrix(rnorm(length(ess_schedule[[1]][[j]])*ncol(path_cur$draws)),
                                                nrow = length(ess_schedule[[1]][[j]]),
                                                ncol = ncol(path_cur$draws))

                        # choose a likelihood threshold
                        threshold <- path_cur$data_log_lik + log(runif(1))

                        # initial proposal, which also defines a bracket
                        theta <- runif(1, 0, 2*pi)
                        lower <- theta - 2*pi; upper <- theta

                        # construct the first proposal
                        if(length(ess_schedule[[1]]) == 1) {
                                copy_mat(draws_prop, cos(theta)*path_cur$draws + sin(theta)*perturbations)
                        } else {
                                copy_2_rows(draws_prop,
                                            cos(theta)*path_cur$draws[ess_schedule[[1]][[j]],] + sin(theta)*perturbations,
                                            ess_schedule[[1]][[j]]-1)
                                copy_2_rows(draws_prop, path_cur$draws[ess_schedule[[2]][[j]],], ess_schedule[[2]][[j]]-1)
                        }

                        # initialize the data log likelihood for the proposed path
                        data_log_lik_prop <- NULL

                        # map the perturbations to an LNA path
                        try({
                                map_draws_2_lna(
                                        pathmat           = pathmat_prop,
                                        draws             = draws_prop,
                                        lna_times         = lna_times,
                                        lna_pars          = lna_parameters,
                                        init_start        = lna_initdist_inds[1],
                                        param_update_inds = param_update_inds,
                                        stoich_matrix     = stoich_matrix,
                                        forcing_inds      = forcing_inds,
                                        forcing_matrix    = forcing_matrix,
                                        svd_sqrt          = svd_sqrt,
                                        svd_d             = svd_d,
                                        svd_U             = svd_U,
                                        svd_V             = svd_V,
                                        lna_pointer       = lna_pointer,
                                        set_pars_pointer  = lna_set_pars_pointer,
                                        step_size         = step_size
                                )

                                census_lna(
                                        path                = pathmat_prop,
                                        census_path         = censusmat,
                                        census_inds         = census_indices,
                                        lna_event_inds      = lna_event_inds,
                                        flow_matrix_lna     = flow_matrix,
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
                                data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                                if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
                        }, silent = TRUE)

                        if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf

                        # continue proposing if not accepted
                        while((data_log_lik_prop < threshold) && !isTRUE(all.equal(lower, upper))) {

                                # increment the number of ESS proposals for the current iteration
                                ess_record[k,j] <- ess_record[k,j] + 1

                                # shrink the bracket
                                if(theta < 0) {
                                        lower <- theta
                                } else {
                                        upper <- theta
                                }

                                # sample a new point
                                theta <- runif(1, lower, upper)

                                # construct the next proposal
                                if(length(ess_schedule[[1]]) == 1) {
                                        copy_mat(draws_prop, cos(theta)*path_cur$draws + sin(theta)*perturbations)
                                } else {
                                        copy_2_rows(draws_prop,
                                                    cos(theta)*path_cur$draws[ess_schedule[[1]][[j]],] +
                                                            sin(theta)*perturbations,
                                                    ess_schedule[[1]][[j]]-1)
                                }

                                # initialize the data log likelihood for the proposed path
                                data_log_lik_prop <- NULL

                                # map the perturbations to an LNA path
                                try({
                                        map_draws_2_lna(
                                                pathmat           = pathmat_prop,
                                                draws             = draws_prop,
                                                lna_times         = lna_times,
                                                lna_pars          = lna_parameters,
                                                init_start        = lna_initdist_inds[1],
                                                param_update_inds = param_update_inds,
                                                stoich_matrix     = stoich_matrix,
                                                forcing_inds      = forcing_inds,
                                                forcing_matrix    = forcing_matrix,
                                                svd_sqrt          = svd_sqrt,
                                                svd_d             = svd_d,
                                                svd_U             = svd_U,
                                                svd_V             = svd_V,
                                                lna_pointer       = lna_pointer,
                                                set_pars_pointer  = lna_set_pars_pointer,
                                                step_size         = step_size
                                        )

                                        census_lna(
                                                path                = pathmat_prop,
                                                census_path         = censusmat,
                                                census_inds         = census_indices,
                                                lna_event_inds      = lna_event_inds,
                                                flow_matrix_lna     = flow_matrix,
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
                                        data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                                        if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
                                }, silent = TRUE)

                                if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
                        }

                        # if the bracket width is not equal to zero, update the draws, path, and data log likelihood
                        if(!isTRUE(all.equal(lower, upper))) {

                                # transfer the new path and residual path into the* sin(theta) path_prop list
                                if(length(ess_schedule[[1]]) == 1) {
                                        copy_mat(path_cur$draws, draws_prop)

                                } else {
                                        copy_2_rows(path_cur$draws,
                                                    draws_prop[ess_schedule[[1]][[j]],],
                                                    ess_schedule[[1]][[j]] - 1)
                                }

                                # copy the LNA path and the data log likelihood
                                copy_mat(path_cur$lna_path, pathmat_prop)
                                copy_vec(path_cur$data_log_lik, data_log_lik_prop)
                        }
                }
        }

        copy_mat(path_cur$ess_record, ess_record)
}