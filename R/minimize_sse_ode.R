#' Obtain estimates of stochastic epidemic model (SEM) parameters on their
#' estimation scales that minimizes the sum of squared errors between the data
#' and the expected path for a SEM approximated by its deterministic ODE limit.
#'
#' NOTE: This should not be considered a reliable method for inference for SEMs.
#' Rather, it is, perhaps, useful as a method for obtaining initial values for
#' other methods.
#'
#' @param stem_object stem_object for which the ODE code and measurement process
#'   have been compiled.
#' @param transformations named list of functions for transforming to and from
#'   the estimation scale, which should be unconstrained (i.e., R^n).
#' @param limits named list with vectors of upper and lower limits for
#'   parameters ON THEIR ESTIMATION SCALE.
#' @param scale either "log" (default) or "linear"
#'
#' @return list containing the parameter estimates on their estimation scales,
#'   Hessian of the data log likelihood evaluated at the maximum likelihood
#'   estimates on their estimation scales, the maximum likelihood ODE path, the
#'   data log-likelihood under the SSE minimizer, the quasi-likelihood variance
#'   inflation factor, and an estimate of the sandwich standard errors for
#'   non-nuisance parameters.
#' @export
minimize_sse_ode <- function(stem_object, transformations = NULL, limits = NULL, scale = "log") {

        # get ode times
        t0        <- stem_object$dynamics$t0
        tmax      <- stem_object$dynamics$tmax
        ode_times <- sort(unique(c(t0, stem_object$dynamics$tcovar[,1], stem_object$measurement_process$data[,1], tmax)))

        # extract the model objects from the stem_object
        popsize                <- stem_object$dynamics$constants["popsize"]
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

        emitmat <- cbind(data[, 1, drop = F],
                         matrix(0.0, nrow = nrow(measproc_indmat), ncol = ncol(measproc_indmat),
                                dimnames = list(NULL, colnames(measproc_indmat))
                         ))

        pathmat_prop <- cbind(ode_times,
                              matrix(0.0, nrow = length(ode_times), ncol = nrow(flow_matrix),
                                     dimnames = list(NULL, c(rownames(flow_matrix)))
                              ))

        # set optimization limits
        lower = pmin(limits$lower, Inf)
        upper = pmin(limits$upper, Inf)

        ### CREATE AUXILLIARY FUNCTIONS ------------------------------------

        # function for getting rid of zeros by the following weighted midpoint rule if scale == "log"
        # let Y_1,...,Y_L be the counts and suppose y_l = 0. Then,
        # imputed_dat = 0.5 * 1/popsize + 0.5 * predicted value from poisson regression glm fit with two points on either side (where available)
        inflate_zeros <- function(dat) {

                which_zeros     <- which(dat[, -1, drop = FALSE] == 0, arr.ind = T)
                which_zeros[,2] <- which_zeros[,2] + 1

                for(r in seq_len(nrow(which_zeros))) {

                        rowinds <- (which_zeros[r,1] - 2):(which_zeros[r,1] + 2)
                        rowinds <- rowinds[rowinds>=1 & rowinds <= nrow(dat)]
                        predmat <- dat[rowinds, c(1,which_zeros[r,2])]

                        pred_vals <- suppressWarnings(exp(predict(glm(predmat[,2] ~ predmat[,1], family = poisson()))))
                        dat[which_zeros[r,1], which_zeros[r,2]] <- 1/popsize +
                                (1 - 1/popsize) * pred_vals[which(predmat[,1] == which_zeros[r,1])]
                }

                return(dat)
        }

        ### inflate zeros if scale == "log"
        if(scale == "log") {
                data <- inflate_zeros(dat = data)
                data[,-1][measproc_indmat] <- log(data[,-1][measproc_indmat])
        }

        # scale transformation functions
        if(is.null(transformations)) {
                transformations <- list(
                        to_estimation_scale   = function(pars) {pars},
                        from_estimation_scale = function(pars) {pars}
                )
        }

        # Compute mean observd incidence
        get_means <- function(pars) {

                pars_nat <- transformations$from_estimation_scale(pars)

                pars2lnapars(ode_pars, pars_nat)
                if(!fixed_inits) init_state <- pars[ode_initdist_inds + 1]

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

                        census_lna(
                                path                = path$incid_path,
                                census_path         = censusmat,
                                census_inds         = census_indices,
                                lna_event_inds      = ode_event_inds,
                                flow_matrix_lna     = t(stoich_matrix),
                                do_prevalence       = do_prevalence,
                                init_state          = init_state,
                                forcing_matrix      = forcing_matrix
                        )

                }, silent = TRUE)

                meanmat <- simulate_r_measure(censusmat       = censusmat,
                                              measproc_indmat = measproc_indmat,
                                              parameters      = pars_nat,
                                              constants       = constants,
                                              tcovar          = tcovar_censmat,
                                              r_measure_ptr   = stem_object$measurement_process$meas_pointers$m_measure_ptr)
                if(scale == "log") meanmat[,-1][measproc_indmat] <- log(meanmat[,-1][measproc_indmat])

                return(meanmat)
        }

        # log-likelihood function
        ode_sse <- function(pars) {
                meanmat <- get_means(pars)
                return(sum((meanmat[,-1][measproc_indmat] - data[,-1][measproc_indmat])^2))
        }

        # compute gradients of the means
        get_grads <- function(pars) {

                grads <- matrix(0.0, nrow = nrow(data) * (ncol(data)-1), ncol = length(pars))
                get_one_mean <- function(pars, j, k) get_means(pars)[j,k+1]

                for(k in seq_len(ncol(measproc_indmat))) {
                        for(j in seq_len(nrow(measproc_indmat))) {
                                grads[(k-1)+j,] <- numDeriv::grad(get_one_mean, pars, j=j,k=k)
                        }
                }

                return(grads)
        }

        # Compute expected variances of observed incidence
        get_vars <- function(pars) {

                pars_nat <- transformations$from_estimation_scale(pars)

                pars2lnapars(ode_pars, pars_nat)
                if(!fixed_inits) init_state <- pars[ode_initdist_inds + 1]

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

                        census_lna(
                                path                = path$incid_path,
                                census_path         = censusmat,
                                census_inds         = census_indices,
                                lna_event_inds      = ode_event_inds,
                                flow_matrix_lna     = t(stoich_matrix),
                                do_prevalence       = do_prevalence,
                                init_state          = init_state,
                                forcing_matrix      = forcing_matrix
                        )

                }, silent = TRUE)

                varmat <- simulate_r_measure(censusmat       = censusmat,
                                             measproc_indmat = measproc_indmat,
                                             parameters      = pars_nat,
                                             constants       = constants,
                                             tcovar          = tcovar_censmat,
                                             r_measure_ptr   = stem_object$measurement_process$meas_pointers$v_measure_ptr)

                if(scale == "log") {
                        varmat[,-1] <- 1 / varmat[,-1]
                }

                return(varmat)
        }

        sandwich <- function(pars, dat, measproc_indmat) {

                datvec <- c(dat[,-1][measproc_indmat])
                grads  <- get_grads(pars)
                means  <- c(get_means(pars)[,-1][measproc_indmat])
                vars   <- c(get_vars(pars)[,-1][measproc_indmat])
                precs  <- diag(vars^-1)
                sighat <- diag((datvec - means)^2)

                # get rid of stupid nuisance parameters in the measurement process
                not_nuisance <- which(apply(grads, 2, function(x) any(x!=0)))
                grads <- grads[,not_nuisance]

                bread <- MASS::ginv(t(grads) %*% precs %*% grads)
                meat  <- t(grads) %*% precs %*% sighat %*% precs %*% grads

                tasty <- bread %*% meat %*% bread
                colnames(tasty) <- rownames(tasty) <- names(pars)[not_nuisance]

                return(tasty)
        }

        ### Get MLEs and Hessian -------------------------------------------
        pars_init <- transformations$to_estimation_scale(stem_object$dynamics$parameters)
        ests <- suppressWarnings(as.numeric(try({
                        optimx::optimx(
                                pars_init,
                                ode_sse,
                                method = c("nlminb"),
                                hessian = FALSE,
                                itnmax = 1e6,
                                upper = upper,
                                lower = lower
                        )}, silent = T)[seq_len(length(pars_init))]))
        hess <- numDeriv::hessian(ode_sse, ests)
        vcov_est <- solve(hess)

        if(any(diag(vcov_est) < 0)) {
                ests <- suppressWarnings(as.numeric(try({
                        optimx::optimx(
                                pars_init,
                                ode_sse,
                                method = c("hjkb"),
                                hessian = FALSE,
                                itnmax = 1e6
                        )}, silent = T)[seq_len(length(pars_init))]))
                hess <- numDeriv::hessian(ode_sse, ests)
                vcov_est <- solve(hess)
        }

        if(any(diag(vcov_est) < 0)) {
                ests <- suppressWarnings(as.numeric(try({
                        optimx::optimx(
                                pars_init,
                                ode_sse,
                                method = c("Nelder-Mead"),
                                hessian = FALSE,
                                itnmax = 1e6
                        )}, silent = T)[seq_len(length(pars_init))]))
                hess <- numDeriv::hessian(ode_sse, ests)
                vcov_est <- solve(hess)
        }

        names(ests) <- colnames(hess) <- rownames(hess) <- names(pars_init)

        # get the ML ODE path and data log likelihood
        ml_pars <- transformations$from_estimation_scale(ests)

        pars2lnapars(ode_pars, ml_pars)
        if(!fixed_inits) init_state <- ml_pars[ode_initdist_inds + 1]

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

                census_lna(
                        path                = path$incid_path,
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

        ### QUASI-LIKELIHOOD -----------------------------------------------
        mean_mat <- get_means(pars = ests)
        var_mat  <- get_vars(pars = ests)

        alpha <- sum((mean_mat[,-1][measproc_indmat] - data[,-1][measproc_indmat])^2 / var_mat[,-1][measproc_indmat]) /
                (sum(measproc_indmat) - length(ml_pars) - 1)

        ### SANDWICH
        sandwich_cov <- sandwich(ests, data, measproc_indmat)

        ### compile estimates and apply delta method
        pars_nat = ml_pars
        pars_est = ests

        nuisance <- which(diag(hess) == 0)

        if(length(nuisance)!= 0) hess <- hess + diag(1e-6, length(pars_est))

        vcov_nat   = mrds::DeltaMethod(par   = pars_est,
                                       fct   = transformations$from_estimation_scale,
                                       vcov  = vcov_est,
                                       delta = 1e-6
                                       )$variance

        rownames(vcov_est) <- rownames(vcov_nat) <- names(pars_est)
        colnames(vcov_est) <- colnames(vcov_nat) <- names(pars_est)

        if(length(nuisance) != 0) {
                pars_est      = pars_est[-nuisance]
                pars_nat      = pars_nat[-nuisance]
                vcov_est      = vcov_est[-nuisance,-nuisance]
                vcov_nat      = vcov_nat[-nuisance,-nuisance]
        }

        ### Return results
        return(
                list(
                        pars_est      = pars_est,
                        pars_nat      = pars_nat,
                        vcov_est      = vcov_est,
                        vcov_nat      = vcov_nat,
                        alpha         = alpha,
                        sandwich_cov  = sandwich_cov,
                        data_log_lik  = data_log_lik,
                        path          = path
                )
        )
}
