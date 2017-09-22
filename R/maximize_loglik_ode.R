#' Obtain estimates of stochastic epidemic model (SEM) parameters that maximize
#' the observed data likelihood for a SEM approximated by its deterministic ODE
#' limit.
#'
#' NOTE: This should not be considered a reliable method for inference for SEMs.
#' Rather, it is perhaps useful as a method for obtaining initial values for
#' other methods.
#'
#' @param stem_object stem_object with
#' @param bounds list containing named vectors of lower and upper bounds for the
#'   model parameters and initial compartment volumes, e.g., \code{list(upper =
#'   c(beta = 1, mu = 1, S_0 = 99, I_0 = 10, R_0 = 10), lower = c(beta = 0, mu =
#'   0, S_0 = 90, I_0 = 1, R_0 = 0))}.
#' @param method optimization method to be passed to \code{optim}.
#'
#' @return list containing the parameter estimates and the inverse of the
#'   Hessian.
#' @export
maximize_loglik_ode <- function(stem_object, bounds, method = "spg") {

        # get ode times
        t0        <- stem_object$dynamics$t0
        tmax      <- stem_object$dynamics$tmax
        ode_times <- sort(unique(c(t0, stem_object$dynamics$tcovar[,1], stem_object$measurement_process$data[,1], tmax)))

        # extract the model objects from the stem_object
        flow_matrix            <- stem_object$dynamics$flow_matrix_ode
        stoich_matrix          <- stem_object$dynamics$stoich_matrix_ode
        ode_pointer            <- stem_object$dynamics$ode_pointers$ode_ptr
        ode_set_pars_pointer   <- stem_object$dynamics$ode_pointers$set_ode_params_ptr
        censusmat              <- stem_object$measurement_process$censusmat
        parameters             <- stem_object$dynamics$parameters
        constants              <- stem_object$dynamics$constants
        n_compartments         <- ncol(flow_matrix)
        n_rates                <- nrow(flow_matrix)
        n_odes                 <- 2*n_rates + n_rates^2
        do_prevalence          <- stem_object$measurement_process$ode_prevalence
        do_incidence           <- stem_object$measurement_process$ode_incidence
        ode_event_inds         <- stem_object$measurement_process$incidence_codes_ode
        n_incidence            <- ifelse(ode_event_inds[1] == -1, 0, length(ode_event_inds))
        state_initializer      <- stem_object$dynamics$state_initializer
        fixed_inits            <- stem_object$dynamics$fixed_inits
        n_strata               <- stem_object$dynamics$n_strata
        ode_initdist_inds      <- stem_object$dynamics$ode_initdist_inds
        initdist_params_cur    <- stem_object$dynamics$initdist_params
        t0                     <- stem_object$dynamics$t0
        t0_fixed               <- stem_object$dynamics$t0_fixed
        step_size              <- stem_object$dynamics$dynamics_args$step_size

        # indices of parameters, constants, and time-varying covariates in the ode_params_* matrices
        ode_param_inds <- seq_along(stem_object$dynamics$param_codes) - 1
        ode_const_inds <- length(ode_param_inds) + seq_along(stem_object$dynamics$const_codes) - 1
        ode_tcovar_inds <- length(ode_param_inds) + length(ode_const_inds) + seq_along(stem_object$dynamics$tcovar_codes) - 1

        # measurement process objects
        data                   <- stem_object$measurement_process$data
        measproc_indmat        <- stem_object$measurement_process$measproc_indmat
        d_meas_pointer         <- stem_object$measurement_process$meas_pointers$d_measure_ptr
        obstimes               <- data[,1]
        obstime_inds           <- stem_object$measurement_process$obstime_inds
        tcovar_censmat         <- stem_object$measurement_process$tcovar_censmat
        census_indices         <- unique(c(0, findInterval(obstimes, ode_times) - 1))

        # generate the matrix of parameters, constants, and time-varying covariates
        ode_pars  <- matrix(0.0,
                            nrow = length(ode_times),
                            ncol = length(stem_object$dynamics$ode_rates$ode_param_codes),
                            dimnames = list(NULL, names(stem_object$dynamics$ode_rates$ode_param_codes)))

        # insert parameters, constants, and time-varying covariates
        parameter_inds <- seq_along(stem_object$dynamics$param_codes)
        constant_inds  <- (length(parameter_inds)+1):(length(parameter_inds) + length(stem_object$dynamics$const_codes))
        tcovar_inds    <- (max(constant_inds)+1):ncol(ode_pars)
        ode_pars[, parameter_inds] <- matrix(stem_object$dynamics$parameters,
                                             nrow = nrow(ode_pars),
                                             ncol = length(parameter_inds), byrow = T)
        ode_pars[, constant_inds]  <- matrix(stem_object$dynamics$constants,
                                             nrow = nrow(ode_pars),
                                             ncol = length(constant_inds), byrow = T)

        # Generate forcing indices and update the ODE parameter matrix with forcings
        forcing_inds <- rep(FALSE, length(ode_times))

        if(!is.null(stem_object$dynamics$dynamics_args$tcovar)) {
                tcovar_rowinds <- match(ode_times, stem_object$dynamics$tcovar[,1])
                ode_pars[tcovar_rowinds, tcovar_inds] <- stem_object$dynamics$tcovar[tcovar_rowinds,-1]

                # zero out forcings if necessary
                if(!is.null(stem_object$dynamics$dynamics_args$forcings)) {

                        # get the forcing indices (supplied in the original tcovar matrix)
                        forcing_inds <- as.logical(match(ode_times,
                                                         stem_object$dynamics$dynamics_args$tcovar[,1],
                                                         nomatch = FALSE))
                        zero_inds    <- !forcing_inds

                        # zero out the tcovar elements corresponding to times with no forcings
                        for(l in seq_along(stem_object$dynamics$dynamics_args$forcings)) {
                                ode_pars[zero_inds, stem_object$dynamics$dynamics_args$forcings[[l]]$tcovar_name] = 0
                        }
                }
        }

        # generate some auxilliary objects
        param_update_inds <- ode_times %in% unique(c(t0, tmax, stem_object$dynamics$tcovar[,1]))

        # retrieve the initial state
        ode_initdist_inds <- stem_object$dynamics$ode_initdist_inds
        init_state        <- ode_pars[1, ode_initdist_inds + 1]

        # generate forcings if they are specified
        forcing_matrix <- matrix(0.0,
                                 nrow = length(ode_times),
                                 ncol = length(stem_object$dynamics$comp_codes),
                                 dimnames = list(NULL, names(stem_object$dynamics$comp_codes)))

        if(!is.null(stem_object$dynamics$dynamics_args$forcings)) {

                for(l in seq_along(stem_object$dynamics$dynamics_args$forcings)) {

                        # insert the flow into the forcing matrix
                        forcing_matrix[forcing_inds, stem_object$dynamics$dynamics_args$forcings[[l]]$from] <-
                                forcing_matrix[forcing_inds, stem_object$dynamics$dynamics_args$forcings[[l]]$from] -
                                stem_object$dynamics$dynamics_args$tcovar[, stem_object$dynamics$dynamics_args$forcings[[l]]$tcovar_name]

                        forcing_matrix[forcing_inds, stem_object$dynamics$dynamics_args$forcings[[l]]$to] <-
                                forcing_matrix[forcing_inds, stem_object$dynamics$dynamics_args$forcings[[l]]$to] +
                                stem_object$dynamics$dynamics_args$tcovar[, stem_object$dynamics$dynamics_args$forcings[[l]]$tcovar_name]

                        # update the adjacency matrix to indicate which rates need to be updated
                        affected_rates <- rep(FALSE, nrow(stem_object$dynamics$tcovar_adjmat))

                        for(n in seq_along(stem_object$dynamics$rates)) {
                                affected_rates[n] = grepl(stem_object$dynamics$dynamics_args$forcings[[l]]$from,
                                                          stem_object$dynamics$rates[[n]]$unparsed) |
                                        grepl(stem_object$dynamics$dynamics_args$forcings[[l]]$to,
                                              stem_object$dynamics$rates[[n]]$unparsed)
                        }

                        stem_object$dynamics$tcovar_adjmat[,stem_object$dynamics$dynamics_args$forcings[[l]]$tcovar_name] <-
                                xor(stem_object$dynamics$tcovar_adjmat[,stem_object$dynamics$dynamics_args$forcings[[l]]$tcovar_name],
                                    affected_rates)
                }
        }

        forcing_matrix <- t(forcing_matrix)

        emitmat <- cbind(stem_object$measurement_process$data[, 1, drop = F],
                         matrix(0.0, nrow = nrow(measproc_indmat), ncol = ncol(measproc_indmat),
                                dimnames = list(NULL, colnames(measproc_indmat))
                         ))

        pathmat_prop <- cbind(ode_times,
                              matrix(0.0, nrow = length(ode_times), ncol = nrow(flow_matrix),
                                     dimnames = list(NULL, c(rownames(flow_matrix)))
                              ))

        # vector of parameters to be optimized over
        optim_pars <- stem_object$dynamics$parameters

        # bounds for model parameters
        lower = bounds$lower
        upper = bounds$upper

        # log-likelihood function
        ode_loglik <- function(pars) {

                pars2lnapars(ode_pars, pars)
                if(!fixed_inits) init_state <- pars[ode_initdist_inds + 1]

                data_log_lik <- -Inf
                try({
                        # integrate the ODEs
                        path <- integrate_odes(ode_times         = ode_times,
                                               ode_pars          = ode_pars,
                                               init_start        = ode_initdist_inds[1],
                                               param_update_inds = param_update_inds,
                                               stoich_matrix     = stoich_matrix,
                                               forcing_inds      = forcing_inds,
                                               forcing_matrix    = forcing_matrix,
                                               step_size         = step_size,
                                               ode_pointer       = ode_pointer,
                                               set_pars_pointer  = ode_set_pars_pointer
                        )

                        path$prev_path <- NULL
                        names(path) <- "ode_path"

                        census_lna(
                                path                = path$ode_path,
                                census_path         = censusmat,
                                census_inds         = census_indices,
                                lna_event_inds      = ode_event_inds,
                                flow_matrix_lna     = t(stoich_matrix),
                                do_prevalence       = do_prevalence,
                                init_state          = init_state,
                                forcing_matrix      = forcing_matrix
                        )

                        # evaluate the density of the incidence counts
                        evaluate_d_measure_LNA(
                                emitmat           = emitmat,
                                obsmat            = stem_object$measurement_process$data,
                                censusmat         = censusmat,
                                measproc_indmat   = measproc_indmat,
                                lna_parameters    = ode_pars,
                                lna_param_inds    = ode_param_inds,
                                lna_const_inds    = ode_const_inds,
                                lna_tcovar_inds   = ode_tcovar_inds,
                                param_update_inds = param_update_inds,
                                census_indices    = census_indices,
                                d_meas_ptr        = d_meas_pointer
                        )

                        # compute the data log likelihood
                        data_log_lik <- sum(emitmat[,-1][measproc_indmat])
                        if(is.nan(data_log_lik)) data_log_lik <- -Inf
                }, silent = TRUE)

                return(-data_log_lik)
        }

        ests <- optim(optim_pars, ode_loglik, method = method, lower = lower, upper = upper, hessian = TRUE)
        hess <- numDeriv::hessian(ode_loglik, ests$par)

        return(list(ests = ests, hessian = hess))
}