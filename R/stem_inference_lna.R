#' Approximate Bayesian inference for a stochastic epidemic model via the linear
#' noise approximation.
#'
#' @param stem_object stochastic epidemic model object with model dynamics, the
#'   measurement process, and a dataset.
#' @param iterations number of MCMC iterations
#' @param prior_density function to compute the prior density for the model
#'   parameters
#' @param kernel list containing the MCMC kernel method, proposal covariance
#'   matrix, and an external pointer for the compiled MCMC kernel function
#' @param t0_kernel output of \link{\code{t0_kernel}}, specifying the RWMH
#'   transition kernel for t0 and the truncated normal distribution prior.
#' @param estimation_scales character vector of estimation scales, each the
#'   result of a
#' @param thin_params thinning interval for posterior parameter samples,
#'   defaults to 1
#' @param thin_latent_proc thinning interval for latent paths, defaults to
#'   ceiling(iterations/100)
#' @param initialization_attempts number of attempts to initialize the latent
#'   path before breaking.
#' @param messages should status messages be generated in an external text file?
#'   If so, the iteration number is printed every 1000th iteration. defaults to
#'   TRUE.
#'
#' @return list with parameter posterior samples and MCMC diagnostics
#' @export
stem_inference_lna <- function(stem_object, iterations, prior_density, kernel, t0_kernel, estimation_scales, thin_params, thin_latent_proc, initialization_attempts = 500, messages, monitor_MCMC) {

        # extract the model objects from the stem_object
        flow_matrix            <- stem_object$dynamics$flow_matrix_lna
        lna_pointer            <- stem_object$dynamics$lna_pointers$lna_pointer
        lna_pointer_ess        <- stem_object$dynamics$lna_pointers$lna_pointer_ess
        lna_set_pars_pointer   <- stem_object$dynamics$lna_pointers$lna_set_pars_ptr
        lna_ess_set_pars_ptr   <- stem_object$dynamics$lna_pointers$lna_ess_set_pars_ptr
        censusmat              <- stem_object$measurement_process$censusmat
        parameters             <- stem_object$dynamics$parameters
        constants              <- stem_object$dynamics$constants
        n_compartments         <- ncol(flow_matrix)
        n_rates                <- nrow(flow_matrix)
        n_odes                 <- 2*n_rates + n_rates^2
        comp_codes             <- stem_object$dynamics$lna_comp_codes
        do_prevalence          <- stem_object$measurement_process$lna_prevalence
        do_incidence           <- stem_object$measurement_process$lna_incidence
        incidence_codes        <- stem_object$measurement_process$incidence_codes_lna
        n_incidence            <- ifelse(incidence_codes[1] == -1, 0, length(incidence_codes))
        census_incidence_codes <- incidence_codes + ncol(censusmat) - (1 + length(incidence_codes))
        names(census_incidence_codes) <- names(incidence_codes)
        state_initializer      <- stem_object$dynamics$state_initializer
        fixed_inits            <- stem_object$dynamics$fixed_inits
        n_strata               <- stem_object$dynamics$n_strata
        lna_initdist_inds      <- stem_object$dynamics$lna_initdist_inds
        initdist_parameters    <- stem_object$dynamics$initdist_params
        t0_fixed               <- stem_object$dynamics$t0_fixed

        # measurement process objects
        data                   <- stem_object$measurement_process$data
        measproc_indmat        <- stem_object$measurement_process$measproc_indmat
        d_meas_pointer         <- stem_object$measurement_process$meas_pointers$d_measure_ptr
        obstimes               <- data[,1]
        obstime_inds           <- stem_object$measurement_process$obstime_inds
        tcovar_censmat         <- stem_object$measurement_process$tcovar_censmat

        # if the initial counts are not fixed, construct the initial distribution prior
        if(!fixed_inits) {
                # function for computing the log prior density for the initial comparment counts
                initdist_prior     <- construct_initdist_prior(state_initializer   = state_initializer,
                                                               n_strata            = n_strata,
                                                               constants           = constants)

                # function for sampling the initial compartment counts (independence sampling)
                initdist_sampler   <- construct_initdist_sampler(state_initializer   = state_initializer,
                                                                 n_strata            = n_strata,
                                                                 constants           = constants)

                # vector for storing the log prior densities for the initial compartment counts
                initdist_log_prior <- double(1 + floor(iterations / thin_params))

                # vector for storing the proposed compartment counts
                initdist_params_prop <- initdist_parameters

        } else {
                initdist_parameters  <- NULL # vector of initial distribution parameters
                initdist_prior       <- NULL # function for computing the prior
                initdist_sampler     <- NULL # function for sampling new values
                initdist_log_prior   <- NULL # vector for storing the log priors
                initdist_params_prop <- NULL # vector for storing the proposed compartment counts
        }

        # get the objects for the MCMC kernel
        kernel_method  <- kernel$method
        kernel_covmat  <- kernel$covmat
        adaptive_rwmh  <- kernel_method == "mvn_adaptive"

        # ensure that the covariance matrix is specified in the same order as the parameters (not including initdist params)
        param_names    <- names(parameters)[!names(parameters) %in% c(names(lna_initdist_inds), "t0")]
        kernel_covmat  <- kernel_covmat[param_names, param_names]

        # set up objects for sampling t0 if it is not fixed
        if(!t0_fixed) {
                t0           <- parameters["t0"]
                t0_prop      <- double(1)
                t0_log_prior <- double(1 + floor(iterations / thin_params))
        } else {
                t0           <- NULL
                t0_prop      <- NULL
                t0_log_prior <- NULL
        }

        # extract the remaining kernel arguments
        if(!adaptive_rwmh) {

                # cache the cholesky decomposition
                kernel_chol_cov <- chol(kernel_covmat)

        } else if(adaptive_rwmh) {

                # unpack the settings for adaptive RWMH
                scale_start      <- kernel$kernel_settings$scale_start
                scale_cooling    <- kernel$kernel_settings$scale_cooling
                max_scaling      <- kernel$kernel_settings$max_scaling
                shape_start      <- kernel$kernel_settings$shape_start
                target           <- kernel$kernel_settings$target
                nugget           <- kernel$kernel_settings$nugget[param_names]
                nugget_weight    <- kernel$kernel_settings$nugget_weight
                opt_scaling      <- 2.38^2 / length(param_names)

                if(length(scale_start) == 0) {
                        scale_start    <- iterations + 2
                        adapt_scale    <- FALSE
                        kernel_scaling <- NULL

                } else {
                        adapt_scale    <- TRUE
                        kernel_scaling <- 1.0
                }

                if(length(shape_start) == 0) {
                        shape_start <- iterations + 2
                        adapt_shape <- FALSE
                        param_means <- double(length(param_names))
                        empirical_covmat <- matrix(0.0, length(param_means), length(param_means))

                } else {
                        adapt_shape <- TRUE
                        param_means <- double(length(param_names))
                        empirical_covmat <- matrix(0.0, length(param_means), length(param_means))
                }

                # initialize objects for monitoring the adaptive MCMC if requested
                if(monitor_MCMC) {
                        # initialize objects for saving the scaling
                        proposal_scalings    <- double(1 + floor(iterations/thin_params))
                        proposal_covariances <- array(0.0, dim = c(nrow(kernel_covmat),
                                                                   ncol(kernel_covmat),
                                                                   1 + floor(iterations/thin_params)))
                        # save initial values
                        if(adapt_scale) proposal_scalings[1]      <- kernel_scaling
                        if(adapt_shape) proposal_covariances[,,1] <- kernel_covmat

                } else {
                        proposal_scalings    <- NULL
                        proposal_covariances <- NULL
                }
        }

        # vectors for storing the model parameters on their natural and estimation scales
        # model_params_nat -- model parameters on their natural scales
        # model_params_est -- model parameters on their estimation scales
        # lna_parameters -- matrix containing all model parameters (including initial count params) for all LNA times
        model_params_nat <- parameters[param_names]
        model_params_est <- double(length(param_names))
        to_estimation_scale(natural_params = model_params_nat, scaled_params = model_params_est, scales = estimation_scales)

        # create analogous vectors for storing the proposed parameter values
        params_prop_nat  <- double(length(param_names))
        params_prop_est  <- double(length(param_names))

        # generate other derived objects
        lna_times         <- sort(unique(c(obstimes, stem_object$dynamics$.dynamics_args$tcovar[,1],
                                           ifelse(t0_fixed, constants["t0"], parameters["t0"]), stem_object$dynamics$tmax)))
        n_times           <- length(lna_times)
        n_census_times    <- length(obstimes)
        param_update_inds <- is.na(match(lna_times, obstimes))
        census_indices    <- findInterval(obstimes, lna_times) - 1

        # matrix for storing the LNA parameters
        lna_parameters    <- matrix(0.0,
                                    nrow = length(lna_times),
                                    ncol = length(stem_object$dynamics$lna_param_codes),
                                    dimnames = list(NULL, names(stem_object$dynamics$lna_param_codes)))
        lna_params_prop   <- matrix(0.0,
                                    nrow = length(lna_times),
                                    ncol = length(stem_object$dynamics$lna_param_codes),
                                    dimnames = list(NULL, names(stem_object$dynamics$lna_param_codes)))

        # insert the lna parameters into the parameter matrix
        pars2lnapars(lna_parameters, c(model_params_nat, t0, initdist_parameters))
        pars2lnapars(lna_params_prop, c(model_params_nat, t0, initdist_parameters))

        # get column indices for constants and time-varying covariates
        const_inds  <- seq_along(stem_object$dynamics$const_codes) + length(stem_object$dynamics$param_codes)
        tcovar_inds <- (max(const_inds)+1):ncol(lna_parameters)

        # insert the constants
        lna_parameters[,const_inds]  <- matrix(stem_object$dynamics$constants,
                                               nrow = nrow(lna_parameters),
                                               ncol = length(const_inds), byrow = T)

        lna_params_prop[,const_inds] <- matrix(stem_object$dynamics$constants,
                                               nrow = nrow(lna_parameters),
                                               ncol = length(const_inds), byrow = T)

        # insert time varying covariates
        if(!is.null(stem_object$dynamics$dynamics_args$tcovar)) {
                tcovar_rowinds                <- findInterval(lna_times, stem_object$dynamics$.dynamics_args$tcovar[,1])
                lna_parameters[, tcovar_inds] <- stem_object$dynamics$.dynamics_args$tcovar[tcovar_rowinds,-1]
                lna_params_prop[,tcovar_inds] <- stem_object$dynamics$.dynamics_args$tcovar[tcovar_rowinds,-1]
        }

        # matrix in which to store the emission probabilities
        emitmat <- cbind(data[, 1, drop = F],
                         matrix(
                                 0.0,
                                 nrow = nrow(measproc_indmat),
                                 ncol = ncol(measproc_indmat),
                                 dimnames = list(NULL, colnames(measproc_indmat))
                         ))

        pathmat <- cbind(lna_times,
                         matrix(
                                 0.0,
                                 nrow = length(lna_times),
                                 ncol = nrow(flow_matrix),
                                 dimnames = list(NULL, c(rownames(flow_matrix)))
                         ))

        # set up MCMC objects
        parameter_samples <-
                matrix(
                        0.0,
                        nrow = 1 + floor(iterations / thin_params),
                        ncol = length(stem_object$dynamics$parameters),
                        dimnames = list(NULL, c(names(model_params_nat), names(t0), names(initdist_parameters)))
                )

        latent_paths      <-
                array(0.0, dim = c(
                        length(lna_times),
                        1 + nrow(flow_matrix),
                        1 + floor(iterations / thin_latent_proc)
                ))

        lna_log_lik       <- double(1 + floor(iterations / thin_params))
        data_log_lik      <- double(1 + floor(iterations / thin_params))
        params_log_prior  <- double(1 + floor(iterations / thin_params))

        # initialize the latent path
        path <- initialize_lna(
                data                    = data,
                lna_parameters          = lna_parameters,
                censusmat               = censusmat,
                emitmat                 = emitmat,
                flow_matrix             = flow_matrix,
                lna_pointer             = lna_pointer,
                lna_set_pars_pointer    = lna_set_pars_pointer,
                lna_times               = lna_times,
                lna_initdist_inds       = lna_initdist_inds,
                param_update_inds       = param_update_inds,
                incidence_codes         = incidence_codes,
                census_incidence_codes  = census_incidence_codes,
                census_indices          = census_indices,
                measproc_indmat         = measproc_indmat,
                obstime_inds            = obstime_inds,
                d_meas_pointer          = d_meas_pointer,
                parameters              = parameters,
                constants               = constants,
                tcovar_censmat          = tcovar_censmat,
                do_prevalence           = do_prevalence,
                do_incidence            = do_incidence,
                initialization_attempts = initialization_attempts
        )

        # compute the log likelihood of the latent LNA path
        path <- lna_density2(
                path              = path,
                lna_times         = lna_times,
                lna_pars          = lna_parameters,
                param_update_inds = param_update_inds,
                flow_matrix       = flow_matrix,
                lna_pointer_ess   = lna_pointer_ess,
                lna_ess_set_pars_ptr = lna_ess_set_pars_ptr
        )

        # set the current LNA log likelihood, data log likelihood, and prior log likelihood
        lna_loglik_cur  <- path$lna_log_lik
        data_loglik_cur <- path$data_log_lik
        params_logprior_cur <- prior_density(c(model_params_nat,t0), model_params_est)

        if(!t0_fixed) {
                t0_logprior_cur <- extraDistr:::cpp_dtnorm(x     = t0,
                                                           mu    = t0_kernel$mean,
                                                           sigma = t0_kernel$sd,
                                                           a     = t0_kernel$lower,
                                                           b     = t0_kernel$upper,
                                                           log_prob = TRUE)
                t0_log_prior[1] <- t0_logprior_cur
        }

        if(!fixed_inits) {
                initdist_logprior_cur <- initdist_prior(initdist_parameters)
                initdist_log_prior[1] <- initdist_logprior_cur
        }

        # save the initial path, data log-likelihood, lna log-likelihood, and prior log-likelihood
        acceptances           <- 0
        path_rec_ind          <- 2 # index for recording the latent paths
        param_rec_ind         <- 2 # index for recording the parameters
        parameter_samples[1,] <- c(model_params_nat, t0, initdist_parameters)
        latent_paths[,,1]     <- path$lna_path
        data_log_lik[1]       <- path$data_log_lik
        lna_log_lik[1]        <- path$lna_log_lik
        params_log_prior[1]   <- params_logprior_cur

        # initialize the status file if status updates are required
        if(messages) {
                status_file <-
                        paste0("LNA_inference_status_",
                               as.numeric(Sys.time()),
                               ".txt")
                cat(
                        "Beginning MCMC",
                        file = status_file,
                        sep = "\n",
                        append = FALSE
                )
        }

        # begin the MCMC
        start.time <- Sys.time()
        for(k in (seq_len(iterations) + 1)) {

                # Print the status if messages are enabled
                if(messages && k%%thin_latent_proc == 0) {
                        # print the iteration
                        cat(paste0("Iteration ", k-1), file = status_file, sep = "\n \n", append = TRUE)
                }

                # Update the path via elliptical slice sampling
                path <- update_lna_path(
                        path_cur                = path,
                        data                    = data,
                        lna_parameters          = lna_parameters,
                        pathmat                 = pathmat,
                        censusmat               = censusmat,
                        emitmat                 = emitmat,
                        flow_matrix             = flow_matrix,
                        lna_pointer_ess         = lna_pointer_ess,
                        lna_ess_set_pars_ptr    = lna_ess_set_pars_ptr,
                        lna_times               = lna_times,
                        lna_initdist_inds       = lna_initdist_inds,
                        param_update_inds       = param_update_inds,
                        incidence_codes         = incidence_codes,
                        census_incidence_codes  = census_incidence_codes,
                        census_indices          = census_indices,
                        measproc_indmat         = measproc_indmat,
                        obstime_inds            = obstime_inds,
                        d_meas_pointer          = d_meas_pointer,
                        parameters              = parameters,
                        constants               = constants,
                        tcovar_censmat          = tcovar_censmat,
                        do_prevalence           = do_prevalence,
                        do_incidence            = do_incidence
                )

                # Sample new parameter values
                if(!adaptive_rwmh) {
                        # propose new parameters via RWMH
                        mvn_rw(params_prop = params_prop_est,
                                                  params_cur = model_params_est,
                                                  cov_chol   = kernel_chol_cov)

                } else {
                        # propose new parameters via adaptive RWMH
                        mvn_adaptive(params_prop      = params_prop_est,
                                     params_cur       = model_params_est,
                                     covmat           = covmat,
                                     empirical_covmat = empirical_covmat,
                                     param_means      = param_means,
                                     scaling          = kernel_scaling,
                                     iteration        = k-1,
                                     acceptances      = acceptances,
                                     nugget           = nugget,
                                     nugget_weight    = nugget_weight,
                                     adapt_scale      = adapt_scale,
                                     adapt_shape      = adapt_shape,
                                     scale_start      = scale_start,
                                     shape_start      = shape_start,
                                     target           = target,
                                     scale_cooling    = scale_cooling,
                                     max_scaling      = max_scaling,
                                     opt_scaling      = opt_scaling)
                }

                # Convert the proposed parameters to their natural scale
                from_estimation_scale(natural_params = params_prop_nat,
                                      scaled_params  = params_prop_est,
                                      scales         = estimation_scales)

                # Compute the log prior for the proposed parameters
                params_logprior_prop <- prior_density(params_prop_nat, params_prop_est)

                # if t0 is not fixed, sample a new value
                if(!t0_fixed) {

                        # sample the new time
                        t0_prop <- extraDistr:::cpp_rtnorm(n     = 1,
                                                           mu    = t0_kernel$mean,
                                                           sigma = t0_kernel$rw_sd,
                                                           a     = t0_kernel$lower,
                                                           b     = t0_kernel$upper)

                        # log prior of the proposal
                        t0_logprior_prop <- extraDistr:::cpp_dtnorm(x        = t0_prop,
                                                                    mu       = t0_kernel$mean,
                                                                    sigma    = t0_kernel$sd,
                                                                    a        = t0_kernel$lower,
                                                                    b        = t0_kernel$upper,
                                                                    log_prob = TRUE)
                        # insert the new time
                        lna_times[1] <- t0_prop
                        path$lna_path[1,1] <- t0_prop
                        path$res_path[1,1] <- t0_prop

                        ### compute the log proposal probabilities for the forward and reverse moves
                        # prob of going from current to new t0
                        t0_cur2new <- extraDistr:::cpp_dtnorm(x         = t0_prop,
                                                              mu        = t0,
                                                              sigma     = t0_kernel$rw_sd,
                                                              a         = t0_kernel$lower,
                                                              b         = t0_kernel$upper,
                                                              log_prob  = TRUE)

                        # prob of going from new to current t0
                        t0_new2cur <- extraDistr:::cpp_dtnorm(x        = t0,
                                                              mu       = t0_prop,
                                                              sigma    = t0_kernel$rw_sd,
                                                              a        = t0_kernel$lower,
                                                              b        = t0_kernel$upper,
                                                              log_prob = TRUE)
                }

                # If the initial compartment counts are not fixed, sample new values and compute the prior
                if(!fixed_inits) {
                        initdist_params_prop    <- initdist_sampler() # sample new initial compartment counts
                        initdist_logprior_prop  <- initdist_prior(initdist_params_prop) # log density of the new counts
                }

                # Insert the proposed parameters into the parameter proposal matrix
                pars2lnapars(lna_params_prop, c(params_prop_nat, t0_prop, initdist_params_prop))

                # Compute the LNA and data log likelihood for the proposed parameters
                path_prop <- lna_density(
                        path              = path,
                        lna_times         = lna_times,
                        lna_pars          = lna_params_prop,
                        param_update_inds = param_update_inds,
                        flow_matrix       = flow_matrix,
                        lna_pointer       = lna_pointer,
                        set_pars_pointer  = lna_set_pars_pointer)

                ## Compute the log posteriors
                # N.B. no need to include the initial distribution log likelihoods
                # since those are updated via an independence sampler so they cancel out
                acceptance_prob <- (path_prop$lna_log_lik + path_prop$data_log_lik + params_logprior_prop) -
                                        (path$lna_log_lik + path$data_log_lik + params_logprior_cur)

                # if t0 is not fixed, need to include the proposal probabilities in the MH ratio
                if(!t0_fixed) acceptance_prob <- acceptance_prob +
                                                        t0_logprior_prop - t0_logprior_cur +
                                                        t0_new2cur - t0_cur2new

                # Accept/Reject via metropolis-hastings (N.B. symmetric parameter proposals)
                if(acceptance_prob >= 0 || acceptance_prob >= log(runif(1))) {

                        ### ACCEPTANCE
                        acceptances         <- acceptances + 1                  # increment acceptances
                        path                <- path_prop                        # update the LNA ODEs
                        lna_parameters      <- lna_params_prop                  # update LNA parameter matrix
                        params_logprior_cur <- params_logprior_prop             # update LNA parameter prior log density
                        copypars(model_params_nat, params_prop_nat)             # update LNA parameters on their natural scales
                        copypars(model_params_est, params_prop_est)             # update LNA parameters on their estimation scales

                        # Update the initial distribution parameters and log prior if not fixed
                        if(!fixed_inits) {
                                initdist_parameters   <- initdist_params_prop   # update initial distribution params
                                initdist_logprior_cur <- initdist_logprior_prop # update initial dist log prior
                        }

                        # update t0 and its log prior if it is not fixed
                        if(!t0_fixed) {
                                t0              <- t0_prop              # update t0
                                t0_logprior_cur <- t0_logprior_prop     # update the log prior for t0
                        }

                } else {
                        ### REJECTION - only need to reset t0 if it is not fixed
                        if(!t0_fixed) {
                                lna_times[1] <- t0
                                path$lna_path[1,1] <- t0
                                path$res_path[1,1] <- t0
                        }
                }

                # Save the latent process if called for in this iteration
                if(k %% thin_latent_proc == 0) {
                        latent_paths[,,path_rec_ind] <- path$lna_path   # save the path
                        path_rec_ind <- path_rec_ind + 1                # increment the path record index
                }

                # Save the parameters if called for in this iteration
                if(k %% thin_params == 0) {

                        # Save the lna log likelihood, data log likelihood, and log priors
                        lna_log_lik[param_rec_ind]      <- path$lna_log_lik
                        data_log_lik[param_rec_ind]     <- path$data_log_lik
                        params_log_prior[param_rec_ind] <- params_logprior_cur

                        if(!t0_fixed)    t0_log_prior[param_rec_ind]       <- t0_logprior_cur
                        if(!fixed_inits) initdist_log_prior[param_rec_ind] <- initdist_logprior_cur

                        # Store the parameter sample
                        parameter_samples[param_rec_ind, ] <- c(model_params_nat, t0, initdist_parameters)

                        # Store the proposal covariance matrix if monitoring is requested
                        if(adaptive_rwmh && monitor_MCMC) {
                                # adaptation_stage = 1, if still using the initial covariance matrix
                                # adaptation_stage = 2, if adapting scaling
                                # adaptation_stage = 3, if adapting shape
                                if(adapt_scale && k >= scale_start && (!adapt_shape || acceptances < shape_start)) {
                                        adaptation_stage <- 2
                                } else if(adapt_shape && acceptances > shape_start) {
                                        adaptation_stage <- 3
                                } else {
                                        adaptation_stage <- 1
                                }

                                # save the scaling and/or proposal covariance
                                if(adaptation_stage == 1) {
                                        if(adapt_scale) proposal_scalings[param_rec_ind] <- 1
                                        if(adapt_shape) proposal_covariances[,,param_rec_ind] <- kernel_covmat

                                        # print the status if asked for
                                        if(monitor_MCMC) cat(paste0("Acceptances = ", acceptances),
                                                             paste0("Acceptance rate = ", acceptances / (k-1)),
                                                             file = status_file, sep = "\n", append = TRUE)

                                } else if(adaptation_stage == 2) {
                                        proposal_scalings[param_rec_ind] <- kernel_scaling
                                        if(adapt_shape) proposal_covariances[,,param_rec_ind] <- kernel_scaling * kernel_covmat

                                        # print the status if asked for
                                        if(monitor_MCMC) cat(paste0("Acceptances = ", acceptances),
                                                             paste0("Acceptance rate = ", acceptances / (k-1)),
                                                             paste0("Scaling factor = ", kernel_scaling),
                                                             file = status_file, sep = "\n", append = TRUE)
                                } else {
                                        if(adapt_scale) proposal_scalings[param_rec_ind] <- opt_scaling
                                        proposal_covariances[,,param_rec_ind] <- empirical_covmat

                                        if(monitor_MCMC) {
                                                cat(paste0("Acceptances = ", acceptances),
                                                    paste0("Acceptance rate = ", acceptances / (k-1)),
                                                    paste0("Scaling factor = ", opt_scaling),
                                                             file = status_file, sep = "\n", append = TRUE)
                                        }
                                }
                        }

                        # Increment the parameter record index
                        param_rec_ind <- param_rec_ind + 1
                }
        }
        end.time <- Sys.time()

        # compile the results
        MCMC_results <- data.frame(lna_log_lik      = lna_log_lik,
                                   data_log_lik     = data_log_lik,
                                   params_log_prior = params_log_prior)

        if(!t0_fixed)    MCMC_results <- cbind(MCMC_results, t0_log_prior = t0_log_prior)
        if(!fixed_inits) MCMC_results <- cbind(MCMC_results, initdist_log_prior = initdist_log_prior)

        MCMC_results <- cbind(MCMC_results, parameter_samples)

        # return the results
        stem_object$dynamics$parameters <- c(model_params_nat, t0, initdist_parameters)
        stem_object$results <- list(time         = difftime(end.time, start.time, units = "hours"),
                                    acceptances  = acceptances,
                                    latent_paths = latent_paths,
                                    MCMC_results = MCMC_results)

        if(monitor_MCMC & adaptive_rwmh) {
                stem_object$results$adaptation_record <- list(proposal_scalings = proposal_scalings,
                                                              proposal_covariances = proposal_covariances)
        }

        return(stem_object)
}