#' Prepare initial distribution list for MCMC
#'
#' @param initializer initializer list
#' @param param_codes vector of parameter codes in the params matrix
#' @param comp_size_vec vector of compartment volumes
#'
#' @return list of objects for MCMC
#' @export
prepare_initdist_objects =
    function(initializer, param_codes, comp_size_vec) {

        # list for initial compartment volume objects
        initdist_objects <- vector("list", length = length(initializer))

        # where in param codes do the initial distributions start
        initdist_start = param_codes[min(grep("_0", names(param_codes)))]

        for(s in seq_along(initdist_objects)) {

            # hyperparameters
            comp_prior <-
                if(!initializer[[s]]$fixed) {
                    if(!is.null(initializer[[s]]$prior)) {
                        if(comp_size_vec[s] != 0) {
                            initializer[[s]]$prior
                        } else {
                            rep(0.0, length(initializer[[s]]$init_states))
                        }
                    } else {
                        if(comp_size_vec[s] != 0) {
                            initializer[[s]]$init_states
                        } else {
                            rep(0.0, length(initializer[[s]]$init_states))
                        }
                    }
                } else {
                    if(comp_size_vec[s] != 0) {
                        initializer[[s]]$init_states
                    } else {
                        rep(0.0, length(initializer[[s]]$init_states))
                    }
                }

            # compartment probabilities
            comp_probs <- if(sum(comp_prior) != 0) {
                comp_prior / sum(comp_prior)
            } else {
                rep(0.0, length(comp_prior))
            }

            # unconstrained moments
            comp_mean  <- comp_size_vec[s] * comp_probs
            comp_cov   <- comp_size_vec[s] * (diag(comp_probs) - comp_probs %*% t(comp_probs))

            if(initializer[[s]]$dist == "dirmultinom")
                comp_cov <- comp_cov * ((comp_size_vec[s] + sum(comp_prior))/(1 + sum(comp_prior)))

            comp_cov_svd <- svd(comp_cov)
            comp_cov_svd$d[length(comp_cov_svd$d)] <- 0
            comp_sqrt_cov <- comp_cov_svd$u %*% diag(sqrt(comp_cov_svd$d))

            # number of compartments
            n_comps_strat = length(initializer[[s]]$init_states)
            initvol_names = paste0(names(initializer[[s]]$init_states),"_0")

            initdist_objects[[s]] <-
                list(init_volumes      = rep(0.0, n_comps_strat),
                     init_volumes_prop = rep(0.0, n_comps_strat),
                     draws_cur         =
                         if(initializer[[s]]$fixed) {
                             rep(0, length(initializer[[s]]$init_states) - 1)
                         } else {
                             rnorm(length(initializer[[s]]$init_states) - 1)
                         },
                     draws_prop        =
                         if(initializer[[s]]$fixed) {
                             rep(0, length(initializer[[s]]$init_states) - 1)
                         } else {
                             rnorm(length(initializer[[s]]$init_states) - 1)
                         },
                     draws_ess         =
                         if(initializer[[s]]$fixed) {
                             rep(0, length(initializer[[s]]$init_states) - 1)
                         } else {
                             rnorm(length(initializer[[s]]$init_states) - 1)
                         },
                     stratum           = initializer[[s]]$strata,
                     fixed             = initializer[[s]]$fixed,
                     comp_size         = comp_size_vec[s],
                     comp_prior        = comp_prior,
                     comp_mean         = comp_mean,
                     comp_sqrt_cov     = comp_sqrt_cov[,-length(comp_mean)])

            # insert the initial compartment volumes
            copy_vec(dest = initdist_objects[[s]]$init_volumes,
                     orig = initializer[[s]]$init_states)

            # indices in the parameter matrix
            initdist_objects[[s]]$param_inds_Cpp =
                param_codes[match(initvol_names, names(param_codes))]
            initdist_objects[[s]]$param_inds_R =
                initdist_objects[[s]]$param_inds_Cpp + 1

            # indices in the MCMC samples matrix
            initdist_objects[[s]]$rec_inds_Cpp =
                initdist_objects[[s]]$param_inds_Cpp - initdist_start
            initdist_objects[[s]]$rec_inds_R =
                initdist_objects[[s]]$rec_inds_Cpp + 1

            # name the initial volume vectors
            names(initdist_objects[[s]]$init_volumes) = initvol_names
            names(initdist_objects[[s]]$init_volumes_prop) = initvol_names
        }

        # map draws to initial volumes
        for(s in seq_along(initdist_objects)) {

            if(!initializer[[s]]$fixed) {

                bad_draws = T

                while(any(bad_draws)) {
                    if (stem_object$dynamics$initializer[[s]]$dist == "rsbln") {
                        orig <- sbln_normal_to_volume(normal_draws = initdist_objects[[s]]$draws_cur,
                                                      stick_means = comp_prior[1:length(initdist_objects[[s]]$draws_cur)],
                                                      stick_sds = comp_prior[(length(initdist_objects[[s]]$draws_cur) + 1):length(comp_prior)],
                                                      stick_size = initdist_objects[[s]]$comp_size)
                    } else {
                        orig <- initdist_objects[[s]]$comp_mean + c(initdist_objects[[s]]$comp_sqrt_cov %*% initdist_objects[[s]]$draws_cur)
                    }

                    # map draws
                    copy_vec(dest = initdist_objects[[s]]$init_volumes,
                             orig = orig)

                    # check boundary conditions
                    bad_draws =
                        any(initdist_objects[[s]]$init_volumes < 0) |
                        any(initdist_objects[[s]]$init_volumes >
                                initdist_objects[[s]]$comp_size)

                    # if outside boundaries, resample
                    if(bad_draws) {
                        initdist_objects[[s]]$draws_cur =
                            rnorm(initdist_objects[[s]]$draws_cur)
                    }
                }
            }
        }

    return(initdist_objects)
}
