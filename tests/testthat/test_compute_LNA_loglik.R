library(stemr)
library(mvtnorm)

context("Compute the likelihood for an LNA path")

test_that("The LNA log-likelihood is computed correctly for a simple model.", {
        skip_on_cran()

        ### Simulate a dataset
        # set up model
        compartments <- c("S","I","R")
        rates <- list(rate("beta * I", "S", "I", incidence = TRUE),
                      rate("mu", "I", "R"))
        state_initializer <- stem_initializer(c(S = 50000, I = 10, R = 50), fixed = FALSE)
        parameters <- c(beta = 0.00001, mu = 1/7, rho = 0.5)
        tcovar <- NULL
        constants <- NULL
        strata <- NULL
        t0 <- 0; tmax <- 52
        timestep <- NULL
        adjacency <- NULL
        messages <- T
        nsim = 1
        census_times = 0:tmax

        # compile dynamics
        dynamics <- stem_dynamics(rates = rates, parameters = parameters, tmax = tmax, state_initializer = state_initializer, compartments=compartments, strata = strata, tcovar = tcovar, messages = TRUE, compile_ode = F, compile_rates = T)

        # compile the measurement process
        emissions <- list(emission("cases", "negbinomial", c("S2I", "S2I * rho"), incidence = TRUE, obstimes = seq(1,tmax,by=1)))

        measurement_process <- stem_measure(emissions = emissions, dynamics = dynamics, messages = T)
        stem_object <- stem(dynamics = dynamics, measurement_process = measurement_process)

        stem_data <- simulate_stem(stem_object = stem_object, method = "gillespie", paths = TRUE, observations = T, tmax = tmax, nsim = 1, census_times = 0:tmax)

        # grab the dataset
        true_path <- stem_data$paths[[1]]
        dat <- stem_data$datasets[[1]]

        # recompile the measurement process
        measurement_process <- stem_measure(emissions = emissions, dynamics = dynamics, data = dat)
        stem_object <- stem(stem_object = stem_object, measurement_process = measurement_process)

        # initialize the inference
        priors <- list(prior(parameter = "beta", distribution = "gamma", hyperpars = c(shape = 1.2e-5, rate = 1), scale = "log"),
                       prior(parameter = "mu", distribution = "gamma", hyperpars = c(shape = 1/8, rate = 1), scale = "log"),
                       prior(parameter = "rho", distribution = "beta", hyperpars = c(shape1 = 95, shape2 = 100), scale = "logit"))

        covmat <- diag(c(1e-6, 1e-6, 5e-2)); rownames(covmat) <- colnames(covmat) <- c("beta", "mu", "rho")

        kernel <- kernel(method = "mvn_rw", covmat = covmat)

        method = "lna"
        iterations = 1000
        thin_params = 1
        thin_latent_proc = 1
        initialization_attempts = 500
        messages = FALSE
        monitor_MCMC = FALSE

        estimation_scales <- sapply(priors, "[[", "scale")
        prior_density     <- construct_priors(priors, param_codes = stem_object$dynamics$param_codes)

        ### Initialize the LNA object
        # extract the model objects from the stem_object
        flow_matrix            <- stem_object$dynamics$flow_matrix_lna
        lna_pointer            <- stem_object$dynamics$lna_pointers$lna_pointer
        lna_pointer_ess        <- stem_object$dynamics$lna_pointers$lna_pointer_ess
        lna_set_pars_pointer   <- stem_object$dynamics$lna_pointers$lna_set_pars_ptr
        lna_ess_set_pars_ptr   <- stem_object$dynamics$lna_pointers$lna_ess_set_pars_ptr
        censusmat              <- stem_object$measurement_process$censusmat
        parameters             <- stem_object$dynamics$parameters
        constants              <- stem_object$dynamics$constants
        initdist_parameters    <- stem_object$dynamics$initdist_params
        lna_initdist_inds      <- stem_object$dynamics$lna_initdist_inds
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

        # measurement process objects
        data                   <- stem_object$measurement_process$data
        measproc_indmat        <- stem_object$measurement_process$measproc_indmat
        d_meas_pointer         <- stem_object$measurement_process$meas_pointers$d_measure_ptr
        obstimes               <- data[,1]
        obstime_inds           <- stem_object$measurement_process$obstime_inds
        tcovar_censmat         <- stem_object$measurement_process$tcovar_censmat

        # get the objects for the MCMC kernel
        kernel_method  <- kernel$method
        kernel_covmat  <- kernel$covmat
        adaptive_rwmh  <- kernel_method == "mvn_adaptive"

        # ensure that the covariance matrix is specified in the same order as the parameters (not including initdist params)
        param_names    <- names(parameters)[!names(parameters) %in% names(lna_initdist_inds)]
        kernel_covmat  <- kernel_covmat[param_names, param_names]

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

                if(is.null(scale_start)) {
                        scale_start    <- iterations + 2
                        adapt_scale    <- FALSE
                        kernel_scaling <- NULL

                } else {
                        adapt_scale    <- TRUE
                        kernel_scaling <- 1
                }

                if(is.null(shape_start)) {
                        shape_start <- iterations + 2
                        adapt_shape <- FALSE
                        param_means <- NULL
                        empirical_covmat <- NULL

                } else {
                        adapt_shape <- TRUE
                        param_means <- parameters[param_names]
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
        model_params_est <- parameters[param_names]
        to_estimation_scale(natural_params = model_params_nat, scaled_params = model_params_est, scales = estimation_scales)

        # create analogous vectors for storing the proposed parameter values
        params_prop_nat  <- parameters[param_names]
        params_prop_est  <- parameters[param_names]

        # generate other derived objects
        lna_times         <- sort(unique(c(obstimes, stem_object$dynamics$.dynamics_args$tcovar[,1],
                                           stem_object$dynamics$t0, stem_object$dynamics$tmax)))
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
        pars2lnapars(lna_parameters, c(model_params_nat, initdist_parameters))

        # get column indices for constants and time-varying covariates
        const_inds             <- seq_along(stem_object$dynamics$const_codes) + length(stem_object$dynamics$param_codes)
        tcovar_inds            <- (max(const_inds)+1):ncol(lna_parameters)

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
                        ncol = length(stem_object$dynamics$parameters)
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
                initdist_prior       <- NULL # function for computing the prior
                initdist_sampler     <- NULL # function for sampling new values
                initdist_log_prior   <- NULL # vector for storing the log priors
                initdist_params_prop <- NULL # vector for storing the proposed compartment counts
        }

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

        ### Compare the LNA density to the density of the LNA computed manually
        lna_loglik_manual <- 0
        for(k in seq_along(lna_times[-1]) + 1) {
                lna_loglik_manual <- lna_loglik_manual + dmvnorm(path$res_path[k,],
                                                                 path$residual[k,],
                                                                 path$diffusion[,,k], log = T)
        }
})