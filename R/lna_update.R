#' Sample a new LNA path via elliptical slice sampling.
#'
#' @inheritParams mvnss_update
#' @param lna_ess_schedule
#'
#' @return list with an updated LNA path along with its stochastic
#'   perturbations, the observed dat log-likelihood, and the lna
#'   log-likelihood, and a record of the number of elliptical slice sampling
#'   proposals
#' @export
lna_update <-
    function(path,
             dat,
             iter,
             parmat,
             lna_ess_schedule,
             lna_ess_control,
             initdist_objects,
             tparam,
             pathmat_prop,
             censusmat,
             draws_prop,
             ess_draws_prop,
             emitmat,
             flow_matrix,
             stoich_matrix,
             census_times,
             forcing_inds,
             forcing_tcov_inds,
             forcings_out,
             forcing_transfers,
             param_vec,
             param_inds,
             const_inds,
             tcovar_inds,
             initdist_inds,
             param_update_inds,
             census_indices,
             event_inds,
             measproc_indmat,
             svd_d,
             svd_U,
             svd_V,
             proc_pointer,
             set_pars_pointer,
             d_meas_pointer,
             do_prevalence,
             joint_initdist_update,
             step_size) {

        # reset ESS steps and angles
        for(s in seq_along(lna_ess_schedule)) {
            reset_vec(lna_ess_schedule[[s]]$steps, 1.0)
            reset_vec(lna_ess_schedule[[s]]$angles, 0.0)
        }

        # perform the elliptical slice sampling updates
        for(k in seq_len(lna_ess_control$n_updates)) {

            # order of strata updates
            ess_order <- sample.int(length(lna_ess_schedule))

            for(j in ess_order) {

                # sample a new set of stochastic perturbations
                ess_draws_prop[lna_ess_schedule[[j]]$ess_inds,] <-
                    rnorm(ess_draws_prop[lna_ess_schedule[[j]]$ess_inds,])

                # choose a likelihood threshold
                threshold <- path$data_log_lik + log(runif(1))

                # initial proposal, which also defines a bracket
                # theta <- runif(1, 0, ess_bracket_width)
                # lower <- theta - ess_bracket_width; upper <- theta
                pos <- runif(1)
                lower <- -lna_ess_schedule[[j]]$bracket_width * pos
                upper <- lower + lna_ess_schedule[[j]]$bracket_width
                theta <- runif(1, lower, upper)

                # initialize the data log likelihood for the proposed path
                data_log_lik_prop <- NULL

                # propose a new initial state
                if(joint_initdist_update) {

                    # indices of strata to sample
                    initdist_codes <- lna_ess_schedule[[j]]$initdist_codes # block in initdist_objects
                    bad_draws      <- vector("logical", length(initdist_codes))

                    # choose an ellipse
                    for(s in initdist_codes) {

                        # if the state is not fixed draw new values
                        if(!initdist_objects[[s]]$fixed) {

                            # draw N(0,1)
                            draw_normals(initdist_objects[[s]]$draws_prop)

                            # compute the linear combination
                            copy_vec(dest = initdist_objects[[s]]$draws_ess,
                                     orig =
                                         cos(theta) * initdist_objects[[s]]$draws_cur +
                                         sin(theta) * initdist_objects[[s]]$draws_prop)

                            if (stem_object$dynamics$initializer[[s]]$dist == "sbln") {
                                orig <- sbln_normal_to_volume(normal_draws = initdist_objects[[s]]$draws_ess,
                                                              stick_means = comp_prior[1:length(initdist_objects[[s]]$draws_cur)],
                                                              stick_sds = comp_prior[(length(initdist_objects[[s]]$draws_cur) + 1):length(comp_prior)],
                                                              stick_size = initdist_objects[[s]]$comp_size)
                            } else {
                                orig <- initdist_objects[[s]]$comp_mean + c(initdist_objects[[s]]$comp_sqrt_cov %*% initdist_objects[[s]]$draws_ess)
                            }


                            # map to volumes
                            copy_vec(dest = initdist_objects[[s]]$init_volumes_prop,
                                     orig = orig)

                            # check boundary conditions
                            bad_draws[s] <-
                                any(initdist_objects[[s]]$init_volumes_prop < 0 |
                                        initdist_objects[[s]]$init_volumes_prop >
                                        initdist_objects[[s]]$comp_size)
                        }
                    }
                }

                if(joint_initdist_update && any(bad_draws)) {
                    data_log_lik_prop <- -Inf

                } else {

                    # copy the new initial compartment counts
                    if(joint_initdist_update) {
                        insert_initdist(parmat = parmat,
                                        initdist_objects = initdist_objects[initdist_codes],
                                        prop = TRUE,
                                        rowind = 0,
                                        mcmc_rec = FALSE)

                        if(!is.null(tparam)) {
                            for(p in seq_along(tparam)) {
                                if(tparam[[p]]$init_dep) {

                                    insert_tparam(
                                        tcovar = parmat,
                                        values =
                                            tparam[[p]]$draws2par(
                                                parameters = parmat[1,],
                                                draws = tparam[[p]]$draws_cur),
                                        col_ind = tparam[[p]]$col_ind,
                                        tpar_inds = tparam[[p]]$tpar_inds_Cpp)
                                }
                            }
                        }
                    }

                    # construct the first proposal
                    copy_2_rows(dest = draws_prop,
                                orig =
                                    cos(theta) * path$draws[lna_ess_schedule[[j]]$ess_inds, ] +
                                    sin(theta) * ess_draws_prop[lna_ess_schedule[[j]]$ess_inds, ],
                                inds = lna_ess_schedule[[j]]$ess_inds - 1)

                    # strata not resampled
                    copy_2_rows(dest = draws_prop,
                                orig = path$draws[lna_ess_schedule[[j]]$complementary_inds,],
                                inds = lna_ess_schedule[[j]]$complementary_inds - 1)

                    try({
                        # map draws onto a latent path
                        map_draws_2_lna(
                            pathmat           = pathmat_prop,
                            draws             = draws_prop,
                            lna_times         = census_times,
                            lna_pars          = parmat,
                            lna_param_vec     = param_vec,
                            lna_param_inds    = param_inds,
                            lna_tcovar_inds   = tcovar_inds,
                            init_start        = initdist_inds[1],
                            param_update_inds = param_update_inds,
                            stoich_matrix     = stoich_matrix,
                            forcing_inds      = forcing_inds,
                            forcing_tcov_inds = forcing_tcov_inds,
                            forcings_out      = forcings_out,
                            forcing_transfers = forcing_transfers,
                            svd_d             = svd_d,
                            svd_U             = svd_U,
                            svd_V             = svd_V,
                            lna_pointer       = proc_pointer,
                            set_pars_pointer  = set_pars_pointer,
                            step_size         = step_size
                        )

                        census_latent_path(
                            path                = pathmat_prop,
                            census_path         = censusmat,
                            census_inds         = census_indices,
                            event_inds          = event_inds,
                            flow_matrix         = flow_matrix,
                            do_prevalence       = do_prevalence,
                            parmat              = parmat,
                            initdist_inds       = initdist_inds,
                            forcing_inds        = forcing_inds,
                            forcing_tcov_inds   = forcing_tcov_inds,
                            forcings_out        = forcings_out,
                            forcing_transfers   = forcing_transfers
                        )

                        # evaluate the density of the incidence counts
                        evaluate_d_measure_LNA(
                            emitmat           = emitmat,
                            obsmat            = dat,
                            censusmat         = censusmat,
                            measproc_indmat   = measproc_indmat,
                            parameters        = parmat,
                            param_inds        = param_inds,
                            const_inds        = const_inds,
                            tcovar_inds       = tcovar_inds,
                            param_update_inds = param_update_inds,
                            census_indices    = census_indices,
                            param_vec         = param_vec,
                            d_meas_ptr        = d_meas_pointer)

                        # compute the data log likelihood
                        data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                        if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
                    }, silent = TRUE)

                    # if proposal failed data_log_lik_prop is -Inf
                    if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
                }

                # continue proposing if not accepted
                while((upper - lower) > sqrt(.Machine$double.eps) &&
                      (data_log_lik_prop < threshold)) {

                    # increment the number of ESS steps
                    increment_elem(lna_ess_schedule[[j]]$steps, k-1)

                    # shrink the bracket
                    if(theta < 0) {
                        lower <- theta
                    } else {
                        upper <- theta
                    }

                    # sample a new point
                    theta <- runif(1, lower, upper)

                    # construct the next initial distribution proposal
                    if(joint_initdist_update) {
                        for(s in initdist_codes) {

                            # if the state is not fixed draw new values
                            if(!initdist_objects[[s]]$fixed) {

                                # compute the linear combination
                                copy_vec(dest = initdist_objects[[s]]$draws_ess,
                                         orig =
                                             cos(theta) * initdist_objects[[s]]$draws_cur +
                                             sin(theta) * initdist_objects[[s]]$draws_prop)

                                if (stem_object$dynamics$initializer[[s]]$dist == "sbln") {
                                    orig <- sbln_normal_to_volume(normal_draws = initdist_objects[[s]]$draws_ess,
                                                                  stick_means = comp_prior[1:length(initdist_objects[[s]]$draws_cur)],
                                                                  stick_sds = comp_prior[(length(initdist_objects[[s]]$draws_cur) + 1):length(comp_prior)],
                                                                  stick_size = initdist_objects[[s]]$comp_size)
                                } else {
                                    orig <- initdist_objects[[s]]$comp_mean + c(initdist_objects[[s]]$comp_sqrt_cov %*% initdist_objects[[s]]$draws_ess)
                                }


                                # map to volumes
                                copy_vec(dest = initdist_objects[[s]]$init_volumes_prop,
                                         orig = orig)

                                # check boundary conditions
                                bad_draws[s] <-
                                    any(initdist_objects[[s]]$init_volumes_prop < 0 |
                                            initdist_objects[[s]]$init_volumes_prop >
                                            initdist_objects[[s]]$comp_size)
                            }
                        }
                    }

                    if(joint_initdist_update && any(bad_draws)) {
                        data_log_lik_prop <- -Inf

                    } else {

                        # copy the new initial compartment counts
                        if(joint_initdist_update) {
                            insert_initdist(parmat = parmat,
                                            initdist_objects = initdist_objects[initdist_codes],
                                            prop = TRUE,
                                            rowind = 0,
                                            mcmc_rec = FALSE)

                            # compute the time-varying parameters if necessary
                            if(!is.null(tparam)) {
                                for(s in seq_along(tparam)) {
                                    if(tparam[[s]]$init_dep) {

                                        # copy into parmat
                                        insert_tparam(
                                            tcovar = parmat,
                                            values =
                                                tparam[[s]]$draws2par(
                                                    parameters = parmat[1,],
                                                    draws = tparam[[s]]$draws_cur),
                                            col_ind = tparam[[s]]$col_ind,
                                            tpar_inds = tparam[[s]]$tpar_inds_Cpp)
                                    }
                                }
                            }
                        }

                        # construct the first proposal
                        copy_2_rows(dest = draws_prop,
                                    orig =
                                        cos(theta) * path$draws[lna_ess_schedule[[j]]$ess_inds, ] +
                                        sin(theta) * ess_draws_prop[lna_ess_schedule[[j]]$ess_inds, ],
                                    inds = lna_ess_schedule[[j]]$ess_inds - 1)

                        # strata not resamples
                        copy_2_rows(dest = draws_prop,
                                    orig = path$draws[lna_ess_schedule[[j]]$complementary_inds,],
                                    inds = lna_ess_schedule[[j]]$complementary_inds - 1)

                        try({
                            # map draws onto a latent path
                            map_draws_2_lna(
                                pathmat           = pathmat_prop,
                                draws             = draws_prop,
                                lna_times         = census_times,
                                lna_pars          = parmat,
                                lna_param_vec     = param_vec,
                                lna_param_inds    = param_inds,
                                lna_tcovar_inds   = tcovar_inds,
                                init_start        = initdist_inds[1],
                                param_update_inds = param_update_inds,
                                stoich_matrix     = stoich_matrix,
                                forcing_inds      = forcing_inds,
                                forcing_tcov_inds = forcing_tcov_inds,
                                forcings_out      = forcings_out,
                                forcing_transfers = forcing_transfers,
                                svd_d             = svd_d,
                                svd_U             = svd_U,
                                svd_V             = svd_V,
                                lna_pointer       = proc_pointer,
                                set_pars_pointer  = set_pars_pointer,
                                step_size         = step_size
                            )

                            census_latent_path(
                                path                = pathmat_prop,
                                census_path         = censusmat,
                                census_inds         = census_indices,
                                event_inds          = event_inds,
                                flow_matrix         = flow_matrix,
                                do_prevalence       = do_prevalence,
                                parmat              = parmat,
                                initdist_inds       = initdist_inds,
                                forcing_inds        = forcing_inds,
                                forcing_tcov_inds   = forcing_tcov_inds,
                                forcings_out        = forcings_out,
                                forcing_transfers   = forcing_transfers
                            )

                            # evaluate the density of the incidence counts
                            evaluate_d_measure_LNA(
                                emitmat           = emitmat,
                                obsmat            = dat,
                                censusmat         = censusmat,
                                measproc_indmat   = measproc_indmat,
                                parameters        = parmat,
                                param_inds        = param_inds,
                                const_inds        = const_inds,
                                tcovar_inds       = tcovar_inds,
                                param_update_inds = param_update_inds,
                                census_indices    = census_indices,
                                param_vec         = param_vec,
                                d_meas_ptr        = d_meas_pointer)

                            # compute the data log likelihood
                            data_log_lik_prop <- sum(emitmat[,-1][measproc_indmat])
                            if(is.nan(data_log_lik_prop)) data_log_lik_prop <- -Inf
                        }, silent = TRUE)

                        # if proposal failed data_log_lik_prop is -Inf
                        if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf
                    }
                }

                # if the bracket width is not equal to zero, update the draws, path, and dat log likelihood
                if((upper - lower) > sqrt(.Machine$double.eps)) {

                    # copy the LNA draws
                    copy_2_rows(dest = path$draws,
                                orig = draws_prop[lna_ess_schedule[[j]]$ess_inds,],
                                inds = lna_ess_schedule[[j]]$ess_inds - 1)

                    # copy the LNA path and the data log likelihood
                    copy_vec(dest = path$data_log_lik, orig = data_log_lik_prop)

                    if(k == lna_ess_control$n_updates & j == length(lna_ess_schedule)) {
                        copy_mat(dest = path$latent_path, orig = pathmat_prop)
                    }

                    # transfer the new initial volumes and draws (volumes already in parameter matrix)
                    if(joint_initdist_update) {
                        for(s in initdist_codes) {
                            if(!initdist_objects[[s]]$fixed) {

                                # copy the N(0,1) draws
                                copy_vec(dest = initdist_objects[[s]]$draws_cur,
                                         orig = initdist_objects[[s]]$draws_ess)

                                # copy the initial compartment volumes
                                copy_vec(dest = initdist_objects[[s]]$init_volumes,
                                         orig = initdist_objects[[s]]$init_volumes_prop)
                            }
                        }

                        # copy time-varying parameters
                        if(!is.null(tparam)) {
                            for(p in seq_along(tparam)) {
                                if(tparam[[p]]$init_dep) {
                                    copy_vec(dest = tparam[[p]]$tpar_cur,
                                             orig = parmat[,tparam[[p]]$col_ind + 1])
                                }
                            }
                        }
                    }

                    # record the final angle
                    insert_elem(dest = lna_ess_schedule[[j]]$angles,
                                elem = theta,
                                ind  = k-1)

                } else {

                    # insert the original compartment counts back into the parameter matrix
                    if(joint_initdist_update) {
                        insert_initdist(parmat = parmat,
                                        initdist_objects = initdist_objects[initdist_codes],
                                        prop = FALSE,
                                        rowind = 0,
                                        mcmc_rec = FALSE)

                        # recover the original time-varying parameter values
                        if(!is.null(tparam)) {
                            for(p in seq_along(tparam)) {
                                if(tparam[[p]]$init_dep) {
                                    vec_2_mat(dest = parmat,
                                              orig = tparam[[p]]$tpar_cur,
                                              ind = tparam[[p]]$col_ind)
                                }
                            }
                        }
                    }
                }
            }
        }

        # update the lna bracket
        if(iter != 0 && iter <= lna_ess_control$bracket_update_iter) {

            for(j in seq_along(lna_ess_schedule)) {

                # angle residual
                copy_vec(
                    dest = lna_ess_schedule[[j]]$angle_resid,
                    orig = mean(lna_ess_schedule[[j]]$angles) -
                        lna_ess_schedule[[j]]$angle_mean)

                # angle variance
                copy_vec(
                    dest = lna_ess_schedule[[j]]$angle_var,
                    orig = lna_ess_schedule[[j]]$angle_resid^2 / iter +
                        lna_ess_schedule[[j]]$angle_var * (iter-1) / iter)

                # angle mean
                copy_vec(
                    dest = lna_ess_schedule[[j]]$angle_mean,
                    orig = mean(lna_ess_schedule[[j]]$angles) / iter +
                        lna_ess_schedule[[j]]$angle_mean * (iter-1) / iter)

                # set the new angle bracket
                if(iter == lna_ess_control$bracket_update_iter) {
                    copy_vec(
                        dest = lna_ess_schedule[[j]]$bracket_width,
                        orig = pmin(lna_ess_control$bracket_scaling *
                                        sqrt(lna_ess_schedule[[j]]$angle_var), 2*pi))
                }
            }
        }
    }
