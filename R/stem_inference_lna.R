#' Approximate Bayesian inference for a stochastic epidemic model via the linear
#' noise approximation.
#'
#' @param stem_object stochastic epidemic model object with model dynamics, the
#'   measurement process, and a dataset.
#' @param iterations number of MCMC iterations
#' @param mcmc_kernel list containing the mcmc_kernel method, proposal
#'   covariance matrix, and an external pointer for the compiled mcmc_kernel
#'   function
#' @param t0_kernel output of \link{\code{t0_kernel}}, specifying the RWMH
#'   transition mcmc_kernel for t0 and the truncated normal distribution prior.
#' @param thin_params thinning interval for posterior parameter samples,
#'   defaults to 1
#' @param thin_latent_proc thinning interval for latent paths, defaults to
#'   ceiling(iterations/100)
#' @param initialization_attempts number of attempts to initialize the latent
#'   path before breaking.
#' @param messages should status messages be generated in an external text file?
#'   If so, the iteration number is printed every 1000th iteration. defaults to
#'   TRUE.
#' @param priors a list of named functions for computing the prior density as
#'   well as transforming parameters to and from their estimation scales. The
#'   functions should have the following names: "prior_density",
#'   "to_estimation_scale", "from_estimation_scale". The prior_density function
#'   must take two vectors as arguments, the model parameters (excluding initial
#'   compartment volumes and t0) on their natural scales, and the model
#'   parameters on their estimation scales. The functions for converting between
#'   parameter scales should take vector of parameters as an argument, returning
#'   a transformed vector (the function call has the form:
#'   \code{transformed_vector <- conversion_function(original_vector)}).
#'   Alternately, the conversion functions may be supplied with two arguments
#'   passed by reference, the first of which is assumed to be the original
#'   vector, and the second which is the transformed vector that will be
#'   modified in place (this case, the function call has the form
#'   \code{conversion_function(original_vector, transformed_vector)}.
#' @param n_ess_updates number of elliptical slice sampling updates per
#'   iteration
#'
#' @return list with parameter posterior samples and MCMC diagnostics
#' @export

stem_inference_lna <- function(stem_object,
                               iterations,
                               priors,
                               mcmc_kernel,
                               t0_kernel,
                               thin_params,
                               thin_latent_proc,
                               initialization_attempts = 500,
                               n_ess_updates = 1,
                               messages) {

        # extract the model objects from the stem_object
        flow_matrix            <- stem_object$dynamics$flow_matrix_lna
        stoich_matrix          <- stem_object$dynamics$stoich_matrix_lna
        lna_pointer            <- stem_object$dynamics$lna_pointers$lna_ptr
        lna_set_pars_pointer   <- stem_object$dynamics$lna_pointers$set_lna_params_ptr
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
        initdist_params_cur    <- stem_object$dynamics$initdist_params
        t0                     <- stem_object$dynamics$t0
        t0_fixed               <- stem_object$dynamics$t0_fixed
        n_ess_updates          <- n_ess_updates

        # measurement process objects
        data                   <- stem_object$measurement_process$data
        measproc_indmat        <- stem_object$measurement_process$measproc_indmat
        d_meas_pointer         <- stem_object$measurement_process$meas_pointers$d_measure_ptr
        obstimes               <- data[,1]
        obstime_inds           <- stem_object$measurement_process$obstime_inds
        tcovar_censmat         <- stem_object$measurement_process$tcovar_censmat

        # construct prior density functions
        prior_density         <- priors$prior_density
        to_estimation_scale   <- priors$to_estimation_scale
        from_estimation_scale <- priors$from_estimation_scale
        in_place_conversion   <- length(formals(to_estimation_scale)) == 2

        # function for converting concentrations to volumes
        if(n_strata == 1) {
                comp_size_vec <- constants["popsize"]
        } else {
                strata_sizes  <- constants[paste0("popsize_", sapply(state_initializer,"[[","strata"))]
                comp_size_vec <- strata_sizes[rep(1:n_strata, sapply(state_initializer, function(x) length(x$codes)))]
        }

        concs2vols <- function(concentrations, size_vec = comp_size_vec) concentrations*size_vec
        vols2concs <- function(volumes, size_vec = comp_size_vec) volumes/size_vec

        # if the initial counts are not fixed, construct the initial distribution prior
        if(!fixed_inits) {
                # function for sampling the initial compartment counts (independence sampling from prior)
                initdist_sampler <- construct_initdist_sampler_lna(state_initializer   = state_initializer,
                                                                   n_strata            = n_strata,
                                                                   constants           = constants)

                # initdist params come in as volumes, convert to concentrations
                init_volumes_cur     <- initdist_params_cur
                init_volumes_prop    <- initdist_params_cur
                initdist_params_cur  <- vols2concs(initdist_params_cur)
                initdist_params_prop <- initdist_params_cur
                names(init_volumes_cur) <- names(init_volumes_prop) <- names(initdist_params_prop) <- names(initdist_params_cur)

        } else {

                init_volumes_cur         <- initdist_params_cur # vector of initial compartment volumes
                initdist_params_cur      <- vols2concs(initdist_params_cur) # vector of initial distribution parameters
                names(init_volumes_cur)  <- names(initdist_params_cur)
                initdist_params_prop     <- NULL # vector for storing the proposed compartment counts
                init_volumes_prop        <- NULL # vector of initial compartment volumes
                initdist_sampler         <- NULL # function for sampling new values
        }

        # grab the names of parameters on their natural and estimation scales
        param_names_nat <- names(parameters)[!names(parameters) %in% c(names(lna_initdist_inds), "t0")]
        param_names_est <- colnames(mcmc_kernel$sigma)
        n_model_params  <- length(param_names_est)

        # vectors for storing the model parameters on their natural and estimation scales
        # model_params_nat -- model parameters on their natural scales
        # model_params_est -- model parameters on their estimation scales
        # lna_params_cur -- matrix containing all model parameters (including initial count params) for all LNA times
        model_params_nat <- parameters[param_names_nat]
        if(in_place_conversion) {
                model_params_est <- double(length(model_params_nat))
                to_estimation_scale(model_params_nat, model_params_est)
        } else {
                model_params_est <- to_estimation_scale(model_params_nat)
        }


        # create analogous vectors for storing the proposed parameter values
        params_prop_nat  <- double(length(param_names_nat)); copy_vec(params_prop_nat, model_params_nat)
        params_prop_est  <- double(length(param_names_nat)); copy_vec(params_prop_est, model_params_est)
        names(params_prop_est) <- names(params_prop_nat) <- param_names_nat

        # generate other derived objects
        lna_times         <- sort(unique(c(obstimes,
                                           stem_object$dynamics$.dynamics_args$tcovar[,1],
                                           stem_object$dynamics$t0,
                                           stem_object$dynamics$tmax)))
        n_times           <- length(lna_times)
        n_census_times    <- length(obstimes)
        param_update_inds <- is.na(match(lna_times, obstimes))
        census_indices    <- findInterval(obstimes, lna_times) - 1

        # set up the MCMC kernel
        if(mcmc_kernel$method == "mvn_rw") {

                acceptances <- 0.0
                sigma_chol  <- chol(mcmc_kernel$sigma)

        } else if(mcmc_kernel$method == "c_rw") {

                acceptances <- rep(0,0, n_model_params)
                kernel_cov  <- diag(mcmc_kernel$sigma)

        } else if(mcmc_kernel$method == "c_rw_adaptive") {

                # MCMC objects
                acceptances      <- rep(0.0, n_model_params)
                adaptations      <- seq_len(iterations)^-mcmc_kernel$kernel_settings$scale_cooling
                proposal_scaling <- rep(1.0, n_model_params)
                std_basis        <- diag(1, n_model_params)

                # empirical mean and covariance of the adaptive kernel
                kernel_mean  <- model_params_est
                kernel_resid <- model_params_est - kernel_mean
                kernel_cov   <- diag(mcmc_kernel$sigma)

                # Adaptation record objects
                adaptation_scale_record <- matrix(1.0,
                                                  nrow = floor(iterations/thin_model_params) + 1,
                                                  ncol = n_model_params)

        } else if(mcmc_kernel$method == "mvn_g_adaptive") {

                # MCMC objects
                acceptances      <- 0.0
                adaptations      <- seq_len(iterations)^-mcmc_kernel$kernel_settings$scale_cooling
                proposal_scaling <- 1.0

                # empirical mean and covariance of the adaptive kernel
                kernel_mean  <- model_params_est
                kernel_resid <- model_params_est - kernel_mean
                kernel_cov   <- mcmc_kernel$sigma

                # Adaptation record objects
                adaptation_scale_record <- matrix(1.0, nrow = floor(iterations/thin_model_params) + 1, ncol = 1)

                adaptation_shape_record <- array(0.0,
                                                  nrow = floor(iterations/thin_model_params) + 1,
                                                  ncol = n_model_params)

        } else if(mcmc_kernel$method == "mvn_c_adaptive") {

                # MCMC objects
                acceptances      <- rep(0.0, n_model_params + 1)
                adaptations      <- seq_len(iterations)^-mcmc_kernel$kernel_settings$scale_cooling
                proposal_scaling <- diag(1, n_model_params)
                std_basis        <- diag(1, n_model_params)

                # empirical mean and covariance of the adaptive kernel
                kernel_mean  <- model_params_est
                kernel_resid <- model_params_est - kernel_mean
                kernel_cov   <- mcmc_kernel$sigma

                # Adaptation record objects
                adaptation_scale_record <- matrix(1.0,
                                                  nrow = floor(iterations/thin_model_params) + 1,
                                                  ncol = n_model_params)

                adaptation_shape_record <- array(0.0,
                                                 nrow = floor(iterations/thin_model_params) + 1,
                                                 ncol = n_model_params)
        }

        # set up objects for sampling t0 if it is not fixed
        if(!t0_fixed) {
                t0_name      <- names(stem_object$dynamics$t0)
                t0_prop      <- double(1)
                t0_log_prior <- double(1 + floor(iterations / thin_params))

                # set the truncation points for the t0 mcmc_kernel if not fixed
                t0_kernel$upper <- min(t0_kernel$upper, min(stem_object$measurement_process$obstimes))
                t0_kernel$lower <- max(t0_kernel$lower, -Inf)
        } else {
                t0           <- NULL
                t0_prop      <- NULL
                t0_log_prior <- NULL
                t0_name      <- NULL
        }

        # matrix for storing the LNA parameters
        lna_params_cur    <- matrix(0.0,
                                    nrow = length(lna_times),
                                    ncol = length(stem_object$dynamics$lna_rates$lna_param_codes),
                                    dimnames = list(NULL, names(stem_object$dynamics$lna_rates$lna_param_codes)))
        lna_params_prop   <- matrix(0.0,
                                    nrow = length(lna_times),
                                    ncol = length(stem_object$dynamics$lna_rates$lna_param_codes),
                                    dimnames = list(NULL, names(stem_object$dynamics$lna_rates$lna_param_codes)))

        # insert the lna parameters into the parameter matrix
        pars2lnapars(lna_params_cur, c(model_params_nat, t0, init_volumes_cur))
        pars2lnapars(lna_params_prop, c(params_prop_nat, t0_prop, init_volumes_prop))

        # get column indices for constants and time-varying covariates
        const_inds  <- seq_along(stem_object$dynamics$const_codes) + length(stem_object$dynamics$param_codes)
        tcovar_inds <- (max(const_inds)+1):ncol(lna_params_cur)

        # insert the constants
        lna_params_cur[,const_inds]  <- matrix(stem_object$dynamics$constants,
                                               nrow = nrow(lna_params_cur),
                                               ncol = length(const_inds), byrow = T)

        lna_params_prop[,const_inds] <- matrix(stem_object$dynamics$constants,
                                               nrow = nrow(lna_params_cur),
                                               ncol = length(const_inds), byrow = T)

        # insert time varying covariates
        if(!is.null(stem_object$dynamics$dynamics_args$tcovar)) {
                tcovar_rowinds                <- findInterval(lna_times, stem_object$dynamics$.dynamics_args$tcovar[,1])
                lna_params_cur[, tcovar_inds] <- stem_object$dynamics$.dynamics_args$tcovar[tcovar_rowinds,-1]
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
        parameter_samples_nat <-
                matrix(
                        0.0,
                        nrow = 1 + floor(iterations / thin_params),
                        ncol = length(stem_object$dynamics$parameters) + !t0_fixed,
                        dimnames = list(NULL, c(names(model_params_nat), t0_name, names(initdist_params_cur)))
                )

        parameter_samples_est <-
                matrix(
                        0.0,
                        nrow = 1 + floor(iterations / thin_params),
                        ncol = length(stem_object$dynamics$parameters) + !t0_fixed,
                        dimnames = list(NULL, c(param_names_est, t0_name, names(initdist_params_cur)))
                )

        lna_paths <- array(0.0, dim = c(n_times, 1 + n_rates, 1 + floor(iterations / thin_latent_proc)))

        data_log_lik      <- double(1 + floor(iterations / thin_params))
        params_log_prior  <- double(1 + floor(iterations / thin_params))
        ess_record        <- matrix(0, nrow = n_ess_updates, ncol = floor(iterations/thin_params))

        # initialize the latent path
        path <- initialize_lna(
                data                    = data,
                lna_parameters          = lna_params_cur,
                censusmat               = censusmat,
                emitmat                 = emitmat,
                stoich_matrix           = stoich_matrix,
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

        # set the current LNA log likelihood, data log likelihood, and prior log likelihood
        params_logprior_cur <- prior_density(model_params_nat, model_params_est)

        if(!t0_fixed) {
                t0_logprior_cur <- extraDistr:::cpp_dtnorm(x        = t0,
                                                           mu       = t0_kernel$mean,
                                                           sigma    = t0_kernel$sd,
                                                           a        = t0_kernel$lower,
                                                           b        = t0_kernel$upper,
                                                           log_prob = TRUE)
                t0_log_prior[1] <- t0_logprior_cur
        }


        # save the initial path, data log-likelihood, lna log-likelihood, and prior log-likelihood
        path_rec_ind          <- 2 # index for recording the latent paths
        param_rec_ind         <- 2 # index for recording the parameters
        parameter_samples_nat[1,] <- c(model_params_nat, t0, initdist_params_cur)
        parameter_samples_est[1,] <- c(model_params_est, t0, init_volumes_cur)
        lna_paths[,,1]        <- path$lna_path
        data_log_lik[1]       <- path$data_log_lik
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
                if((messages) && k%%thin_latent_proc == 0) {
                        # print the iteration
                        cat(paste0("Iteration ", k-1), file = status_file, sep = "\n \n", append = TRUE)
                }

                # Update the path via elliptical slice sampling
                path <- update_lna_path(
                        path_cur                = path,
                        data                    = data,
                        lna_parameters          = lna_params_cur,
                        pathmat                 = pathmat,
                        censusmat               = censusmat,
                        emitmat                 = emitmat,
                        flow_matrix             = flow_matrix,
                        stoich_matrix           = stoich_matrix,
                        lna_pointer             = lna_pointer,
                        lna_set_pars_ptr        = lna_set_pars_pointer,
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
                        n_ess_updates           = n_ess_updates
                )

                # compute the current log posterior
                logpost_cur <- path$data_log_lik + params_logprior_cur

                # Sample new parameter values
                if(mcmc_kernel$method == "c_rw") {

                        # order in which the parameters should be updated
                        component_order <- sample.int(n_model_params, replace = FALSE)

                        for(s in component_order) {

                                # sample the model parameter - s-1 gives C++ style indexing
                                c_rw(params_prop_est, model_params_est, s-1, kernel_cov)

                                # Convert the proposed parameters to their natural scale
                                if(in_place_conversion) {
                                        from_estimation_scale(params_prop_est, params_prop_nat)
                                } else {
                                        params_prop_nat <- from_estimation_scale(params_prop_est)
                                }

                                # Compute the log prior for the proposed parameters
                                params_logprior_prop <- prior_density(params_prop_nat, params_prop_est)

                                # Insert the proposed parameters into the parameter proposal matrix
                                pars2lnapars(lna_params_prop, c(params_prop_nat, t0_prop, init_volumes_prop))

                                # set the data log likelihood for the proposal to NULL
                                data_log_lik_prop <- NULL

                                try({
                                        map_draws_2_lna(pathmat           = path$lna_path_prop,
                                                        draws             = path$draws,
                                                        lna_times         = lna_times,
                                                        lna_pars          = lna_params_prop,
                                                        init_start        = lna_initdist_inds[1],
                                                        param_update_inds = param_update_inds,
                                                        stoich_matrix     = stoich_matrix,
                                                        lna_pointer       = lna_pointer,
                                                        set_pars_pointer  = lna_set_pars_pointer
                                        )

                                        # census the LNA
                                        census_lna(
                                                path                = path$lna_path_prop,
                                                census_path         = censusmat,
                                                census_inds         = census_indices,
                                                flow_matrix_lna     = flow_matrix,
                                                do_prevalence       = do_prevalence,
                                                init_state          = lna_params_prop[1, lna_initdist_inds + 1],
                                                incidence_codes_lna = incidence_codes
                                        )

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
                                }, silent = TRUE)

                                if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf

                                # compute the log posterior for the proposed parameters
                                logpost_prop <- data_log_lik_prop + params_logprior_prop

                                ## Compute the acceptance probability
                                acceptance_prob <- logpost_prop - logpost_cur

                                # Accept/Reject via metropolis-hastings
                                if(acceptance_prob >= 0 || acceptance_prob >= log(runif(1))) {

                                        ### ACCEPTANCE
                                        acceptances[s]      <- acceptances[s] + 1    # increment acceptances

                                        path$data_log_lik   <- data_log_lik_prop     # save the data log likelihood
                                        params_logprior_cur <- params_logprior_prop  # update LNA parameter prior log density
                                        logpost_cur         <- logpost_prop          # save the log posterior

                                        copy_col(lna_params_cur, lna_params_prop, s-1) # update the model parameter
                                        copy_elem(model_params_nat, params_prop_nat, s-1)  # update LNA parameters on their natural scales
                                        copy_elem(model_params_est, params_prop_est, s-1)  # update LNA parameters on their estimation scales
                                }
                        }

                } else if(mcmc_kernel$method == "mvn_rw") {

                        # propose new parameters
                        mvn_rw(params_prop_est, model_params_est, sigma_chol)

                        # Convert the proposed parameters to their natural scale
                        if(in_place_conversion) {
                                from_estimation_scale(params_prop_est, params_prop_nat)
                        } else {
                                params_prop_nat <- from_estimation_scale(params_prop_est)
                        }

                        # Compute the log prior for the proposed parameters
                        params_logprior_prop <- prior_density(params_prop_nat, params_prop_est)

                        # Insert the proposed parameters into the parameter proposal matrix
                        pars2lnapars(lna_params_prop, c(params_prop_nat, t0_prop, init_volumes_prop))

                        # set the data log likelihood for the proposal to NULL
                        data_log_lik_prop <- NULL

                        try({
                                map_draws_2_lna(pathmat           = path$lna_path_prop,
                                                draws             = path$draws,
                                                lna_times         = lna_times,
                                                lna_pars          = lna_params_prop,
                                                init_start        = lna_initdist_inds[1],
                                                param_update_inds = param_update_inds,
                                                stoich_matrix     = stoich_matrix,
                                                lna_pointer       = lna_pointer,
                                                set_pars_pointer  = lna_set_pars_pointer
                                )

                                # census the LNA
                                census_lna(
                                        path                = path$lna_path_prop,
                                        census_path         = censusmat,
                                        census_inds         = census_indices,
                                        flow_matrix_lna     = flow_matrix,
                                        do_prevalence       = do_prevalence,
                                        init_state          = lna_params_prop[1, lna_initdist_inds + 1],
                                        incidence_codes_lna = incidence_codes
                                )

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
                        }, silent = TRUE)

                        if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf

                        # compute the log posterior for the proposed parameters
                        logpost_prop <- data_log_lik_prop + params_logprior_prop

                        ## Compute the acceptance probability
                        acceptance_prob <- logpost_prop - logpost_cur

                        # Accept/Reject via metropolis-hastings
                        if(acceptance_prob >= 0 || acceptance_prob >= log(runif(1))) {

                                ### ACCEPTANCE
                                acceptances         <- acceptances + 1      # increment acceptances

                                path$data_log_lik   <- data_log_lik_prop    # save the data log likelihood
                                params_logprior_cur <- params_logprior_prop # update LNA parameter prior log density
                                logpost_cur         <- logpost_prop         # save the log posterior

                                copy_mat(lna_params_cur, lna_params_prop)   # update the model parameter
                                copy_vec(model_params_nat, params_prop_nat) # update LNA parameters on their natural scales
                                copy_vec(model_params_est, params_prop_est) # update LNA parameters on their estimation scales
                        }

                } else if(mcmc_kernel$method == "c_rw_adaptive") {

                        # order in which the parameters should be updated
                        component_order <- sample.int(n_model_params, replace = FALSE)

                        for(s in component_order) {

                                # sample the model parameter - s-1 gives C++ style indexing
                                c_rw_adaptive(params_prop_est,
                                              model_params_est,
                                              s - 1,
                                              kernel_cov,
                                              proposal_scaling,
                                              adaptations,
                                              nugget)

                                # Convert the proposed parameters to their natural scale
                                if(in_place_conversion) {
                                        from_estimation_scale(params_prop_est, params_prop_nat)
                                } else {
                                        params_prop_nat <- from_estimation_scale(params_prop_est)
                                }

                                # Compute the log prior for the proposed parameters
                                params_logprior_prop <- prior_density(params_prop_nat, params_prop_est)

                                # Insert the proposed parameters into the parameter proposal matrix
                                pars2lnapars(lna_params_prop, c(params_prop_nat, t0_prop, init_volumes_prop))

                                # set the data log likelihood for the proposal to NULL
                                data_log_lik_prop <- NULL

                                try({
                                        map_draws_2_lna(pathmat           = path$lna_path_prop,
                                                        draws             = path$draws,
                                                        lna_times         = lna_times,
                                                        lna_pars          = lna_params_prop,
                                                        init_start        = lna_initdist_inds[1],
                                                        param_update_inds = param_update_inds,
                                                        stoich_matrix     = stoich_matrix,
                                                        lna_pointer       = lna_pointer,
                                                        set_pars_pointer  = lna_set_pars_pointer
                                        )

                                        # census the LNA
                                        census_lna(
                                                path                = path$lna_path_prop,
                                                census_path         = censusmat,
                                                census_inds         = census_indices,
                                                flow_matrix_lna     = flow_matrix,
                                                do_prevalence       = do_prevalence,
                                                init_state          = lna_params_prop[1, lna_initdist_inds + 1],
                                                incidence_codes_lna = incidence_codes
                                        )

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
                                }, silent = TRUE)

                                if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf

                                # compute the log posterior for the proposed parameters
                                logpost_prop <- data_log_lik_prop + params_logprior_prop

                                ## Compute the acceptance probability
                                acceptance_prob <- logpost_prop - logpost_cur

                                # Accept/Reject via metropolis-hastings
                                if(acceptance_prob >= 0 || acceptance_prob >= log(runif(1))) {

                                        ### ACCEPTANCE
                                        acceptances[s]      <- acceptances[s] + 1    # increment acceptances

                                        path$data_log_lik   <- data_log_lik_prop     # save the data log likelihood
                                        params_logprior_cur <- params_logprior_prop  # update LNA parameter prior log density
                                        logpost_cur         <- logpost_prop          # save the log posterior

                                        copy_col(lna_params_cur, lna_params_prop, s-1) # update the model parameter
                                        copy_elem(model_params_nat, params_prop_nat, s-1)  # update LNA parameters on their natural scales
                                        copy_elem(model_params_est, params_prop_est, s-1)  # update LNA parameters on their estimation scales
                                }
                        }

                } else if(mcmc_kernel$method == "mvn_c_adaptive") {

                } else if(mcmc_kernel$method == "mvn_g_adaptive") {

                }

                # Convert the proposed parameters to their natural scale
                params_prop_nat <- from_estimation_scale(params_prop_est)

                # Compute the log prior for the proposed parameters
                params_logprior_prop <- prior_density(params_prop_nat, params_prop_est)

                # if t0 is not fixed, sample a new value
                if(!t0_fixed) {

                        # sample the new time from proposal centered at t0
                        t0_prop <- extraDistr:::cpp_rtnorm(n     = 1,
                                                           mu    = t0,
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
                        initdist_params_prop    <- initdist_sampler() # sample new initial compartment concentrations
                        init_volumes_prop       <- concs2vols(initdist_params_prop) # convert to volumes
                }

                # Insert the proposed parameters into the parameter proposal matrix
                pars2lnapars(lna_params_prop, c(params_prop_nat, t0_prop, init_volumes_prop))

                # map the perturbations to an LNA path and compute the data log likelihood for the proposed parameters
                # it is possible that, with the non-centered parameterization the draws result in a degenerate mapping

                # set the data log likelihood for the proposal to NULL
                data_log_lik_prop <- NULL

                try({
                        map_draws_2_lna(pathmat           = path$lna_path_prop,
                                        draws             = path$draws,
                                        lna_times         = lna_times,
                                        lna_pars          = lna_params_prop,
                                        init_start        = lna_initdist_inds[1],
                                        param_update_inds = param_update_inds,
                                        stoich_matrix     = stoich_matrix,
                                        lna_pointer       = lna_pointer,
                                        set_pars_pointer  = lna_set_pars_pointer
                        )

                        # census the LNA
                        census_lna(
                                path                = path$lna_path_prop,
                                census_path         = censusmat,
                                census_inds         = census_indices,
                                flow_matrix_lna     = flow_matrix,
                                do_prevalence       = do_prevalence,
                                init_state          = lna_params_prop[1, lna_initdist_inds + 1],
                                incidence_codes_lna = incidence_codes
                        )

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
                }, silent = TRUE)

                if(is.null(data_log_lik_prop)) data_log_lik_prop <- -Inf

                # compute the log posterior for the proposed parameters
                logpost_prop <- data_log_lik_prop + params_logprior_prop

                ## Compute the log posteriors
                # N.B. no need to include the initial distribution log likelihoods
                # since those are updated via an independence sampler so they cancel out
                acceptance_prob <- logpost_prop - logpost_cur

                # if t0 is not fixed, need to include the proposal probabilities in the MH ratio
                if(!t0_fixed) acceptance_prob <- acceptance_prob +
                                                        t0_logprior_prop - t0_logprior_cur +
                                                        t0_new2cur - t0_cur2new

                # Accept/Reject via metropolis-hastings
                if(acceptance_prob >= 0 || acceptance_prob >= log(runif(1))) {

                        ### ACCEPTANCE
                        acceptances         <- acceptances + 1                  # increment acceptances
                        path$data_log_lik   <- data_log_lik_prop                # save the data log likelihood
                        params_logprior_cur <- params_logprior_prop             # update LNA parameter prior log density
                        copy_mat(lna_params_cur, lna_params_prop)               # update LNA parameter matrix
                        copy_vec(model_params_nat, params_prop_nat)             # update LNA parameters on their natural scales
                        copy_vec(model_params_est, params_prop_est)             # update LNA parameters on their estimation scales

                        # Update the initial distribution parameters and log prior if not fixed
                        if(!fixed_inits) {
                                initdist_params_cur <- initdist_params_prop   # update initial distribution params
                                init_volumes_cur    <- init_volumes_prop      # update initial volumes
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
                                path$lna_path_prop[1,1] <- t0
                        }
                }

                # Save the latent process if called for in this iteration
                if(k %% thin_latent_proc == 0) {
                        ess_record[,param_rec_ind-1] <- path$ess_record  # save the ESS record
                        lna_paths[,,path_rec_ind]    <- path$lna_path    # save the path
                        path_rec_ind                 <- path_rec_ind + 1 # increment the path record index
                }

                # Save the parameters if called for in this iteration
                if(k %% thin_params == 0) {

                        # Save the lna log likelihood, data log likelihood, and log priors
                        data_log_lik[param_rec_ind]     <- path$data_log_lik
                        params_log_prior[param_rec_ind] <- params_logprior_cur

                        if(!t0_fixed)    t0_log_prior[param_rec_ind]       <- t0_logprior_cur

                        # Store the parameter sample
                        parameter_samples_nat[param_rec_ind, ] <- c(model_params_nat, t0, initdist_params_cur)
                        parameter_samples_est[param_rec_ind, ] <- c(model_params_est, t0, init_volumes_cur)

                        # Store the proposal covariance matrix if monitoring is requested
                        if(adaptive_rwmh) {
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
                                        if(messages) cat(paste0("Acceptances = ", acceptances),
                                                             paste0("Acceptance rate = ", acceptances / (k-1)),
                                                             file = status_file, sep = "\n", append = TRUE)

                                } else if(adaptation_stage == 2) {
                                        proposal_scalings[param_rec_ind] <- kernel_scaling
                                        if(adapt_shape) proposal_covariances[,,param_rec_ind] <- kernel_scaling * kernel_covmat

                                        # print the status if asked for
                                        if(messages) cat(paste0("Acceptances = ", acceptances),
                                                             paste0("Acceptance rate = ", acceptances / (k-1)),
                                                             paste0("Scaling factor = ", kernel_scaling),
                                                             file = status_file, sep = "\n", append = TRUE)
                                } else {
                                        if(adapt_scale) proposal_scalings[param_rec_ind] <- opt_scaling
                                        proposal_covariances[,,param_rec_ind] <- empirical_covmat

                                        if(messages) {
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
        MCMC_results <- data.frame(
                # lna_log_lik      = lna_log_lik, # NOT NEEDED
                data_log_lik     = data_log_lik,
                params_log_prior = params_log_prior
                )

        if(!t0_fixed)    MCMC_results <- cbind(MCMC_results, t0_log_prior = t0_log_prior)
        # if(!fixed_inits) MCMC_results <- cbind(MCMC_results, initdist_log_prior = initdist_log_prior)

        MCMC_results <- cbind(MCMC_results, parameter_samples_nat, parameter_samples_est)

        # coerce the ESS record to a vector if there is only one ESS update per iteration
        if(n_ess_updates == 1) ess_record <- as.vector(ess_record)

        # return the results
        stem_object$dynamics$parameters <- c(model_params_nat, t0, init_volumes_cur)
        names(stem_object$dynamics$parameters) <- names(stem_object$dynamics$param_codes)
        stem_object$results <- list(time         = difftime(end.time, start.time, units = "hours"),
                                    acceptances  = acceptances,
                                    lna_paths = lna_paths,
                                    MCMC_results = MCMC_results,
                                    ess_record   = ess_record)

        # save the settings
        stem_object$stem_settings <- list(iterations       = iterations,
                                          thin_params      = thin_params,
                                          thin_latent_proc = thin_latent_proc,
                                          n_ess_updates    = n_ess_updates,
                                          priors           = priors,
                                          prior_density    = prior_density,
                                          mcmc_kernel      = mcmc_kernel,
                                          t0_kernel        = t0_kernel,
                                          initdist_sampler = initdist_sampler)


        if(adaptive_rwmh) {
                if(adapt_scale && k >= scale_start && (!adapt_shape || acceptances < shape_start)) {
                        stem_object$results$scaling <- kernel_scaling
                } else if(adapt_shape && acceptances > shape_start) {
                        stem_object$results$scaling <- opt_scaling
                } else {
                        stem_object$results$scaling <- 1
                }

                if(adapt_shape) stem_object$results$empirical_covmat <- empirical_covmat

                stem_object$results$adaptation_record <- list(proposal_scalings = proposal_scalings,
                                                              proposal_covariances = proposal_covariances,
                                                              nugget = nugget,
                                                              nugget_weight = nugget_weight)
        }

        return(stem_object)
}