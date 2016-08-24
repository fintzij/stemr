#' Simulations from a stochastic epidemic model.
#'
#' @param stem Stochastic epidemic model object. The dynamics slot must be
#'   supplied at a minimum. If simulations from the measurement process are
#'   desired, the measurement process slot must also be supplied.
#' @param nsim number of realizations to simulate
#' @param paths Should population-level paths be returned?
#' @param observations Should simulated observations be returned? Requires that
#'   a measurement process be defined in the stem object.
#' @param subject_paths Should population-level paths be mapped to subject-level
#'   paths and returned (e.g. for the purpose of initializing a subject-level
#'   collection of disease histories)? Only available for exact simulation via
#'   the Gillespie direct method.
#' @param method either "gillespie" if simulating via Gillespie's direct method,
#'   or "lna" if simulating paths via the Linear Noise Approximation.
#' @param lna_restart Either a logical indicating whether the non-restarting
#'   version of the LNA (default) should be used, or whether the LNA should be
#'   restarted at census times. Alternately, a vector of restart times could be
#'   supplied. NOTE: THE RESTARTING LNA IMPLEMENTATION HAS NUMERICAL PROBLEMS
#'   THAT HAVE NOT YET BEEN FIXED.
#' @param t0 the time at which the system is initialized. If not supplied, it is
#'   taken to be the first observation time if that is available, or otherwise
#'   0. If supplied, it must be no later than the first observation time.
#' @param tmax the time at which simulation of the system is terminated. If not
#'   supplied, defaults to the last observation time if not supplied.
#' @param timestep time discretization for smooth time-varying functions,
#'   specifically seasonal terms, to be used in simulation when
#'   method="gillespie". If not supplied, \code{timestep} will default to
#'   \code{(tmax - t0)/50}.
#' @param census_times vector of times at which compartment counts should be
#'   recorded. Required for \code{method = "lna"}. If supplied and \code{method
#'   = "gillespie"}, the compartment counts at census times are returned rather
#'   than the full paths.
#' @param paths_as_array if TRUE, the simulated paths will be returned as an
#'   array. Available only if \code{census_times} is specified.
#' @param datasets_as_array if TRUE, the simulated datasets will be returned as
#'   an array.
#' @param messages should a message be printed when parsing the rates?
#'
#' @return Returns a list with the simulated paths, subject-level paths, and/or
#'   datasets. If \code{paths = FALSE} and \code{observations = FALSE}, or if
#'   \code{paths = TRUE} and \code{observations = TRUE}, a list \code{nsim stem}
#'   paths and datasets, each returned either as a list or array, is returned.
#'   If \code{subject_paths = TRUE} and \code{method = "gillespie"}, a list of
#'   subject-level paths is also returned.
#'
#'   If \code{paths = TRUE} and \code{observations = FALSE}, a list or array of
#'   simulated population-level paths is returned. If \code{method =
#'   "gillespie"} and \code{subject_paths = TRUE}, an additional list is
#'   returned, containing lists of subject-level mappings. Each sublist contains
#'   a vector of subject states at t0, along with a mapped path.
#'
#'   If \code{paths = FALSE} and \code{observations = TRUE}, a list or array of
#'   simulated datasets is returned.
#' @export
simulate_stem <- function(stem_object, nsim = 1, paths = FALSE, observations = FALSE, subject_paths = FALSE, method = "gillespie", lna_restart = FALSE, t0 = NULL, tmax = NULL, timestep = NULL, census_times = NULL, paths_as_array = FALSE, datasets_as_array = FALSE, messages = TRUE) {

        # ensure that the method is correctly specified
        if(!method %in% c("gillespie", "lna")) {
                stop("The simulation method must either be 'gillespie' or 'lna'.")
        }

        # if lna, subject paths are not available
        if(subject_paths && (method == "lna")) {
                warning("Subject-paths are only available when simulating paths when method='gillespie'.")
        }

        # check that the stem_dynamics are supplied
        if(is.null(stem_object$dynamics)) {
                stop("The stochastic epidemic model dynamics must be specified.")
        }

        # check that a measurement process is supplied if one is required
        if(is.null(stem_object$measurement_process) & ((!paths & !observations) | observations)) {
                stop("In order to simulate a dataset, a measurement process must be supplied in the stem_object.")
        }

        if(paths == observations) {
                paths <- TRUE
                observations <- TRUE
        }


        # build the time varying covariate matrix (includes, at a minimum, the endpoints of the simulation interval)
        # if timestep is null, there are no time-varying covariates
        if(method == "gillespie") {

                # if any of t0, tmax, or a timestep was supplied, check if they differ from the parameters supplied in the stem_object$dynamics.
                # if they differ, reconstruct the tcovar matrix and associated objects
                rebuild_tcovar <- (!is.null(t0) && t0 != stem_object$dynamics$tcovar[1,1]) ||
                        (!is.null(tmax) && tmax != stem_object$dynamics$tcovar[nrow(stem_object$dynamics$tcovar), 1]) || !is.null(timestep)

                if(rebuild_tcovar) {
                        # get t0 or tmax if it wasn't supplied
                        if(is.null(t0)) t0 <- stem_object$dynamics$tcovar[1,1]
                        if(is.null(tmax)) tmax <- stem_object$dynamics$tcovar[nrow(stem_object$dynamics$tcovar), 1]

                        # replace t0 and tmax in the time-varying covariate matrix
                        if(is.null(timestep) && stem_object$dynamics$timevarying) {
                                timestep <- (tmax - t0)/50

                        } else if(is.null(timestep) && !stem_object$dynamics$timevarying){
                                timestep <- tmax - t0
                        }

                        # rebuild the time-varying covariate matrix so that it contains the census intervals
                        stem_object$dynamics$tcovar              <- build_tcovar_matrix(tcovar = stem_object$dynamics$.dynamics_args$tcovar,
                                                                                        timestep = timestep, t0 = t0, tmax = tmax)
                        stem_object$dynamics$tcovar_codes        <- seq_len(ncol(stem_object$dynamics$tcovar) - 1)
                        names(stem_object$dynamics$tcovar_codes) <- colnames(stem_object$dynamics$tcovar)[2:ncol(stem_object$dynamics$tcovar)]
                        stem_object$dynamics$n_tcovar            <- ncol(stem_object$dynamics$tcovar) - 1
                        stem_object$dynamics$tcovar_changemat    <- build_tcovar_changemat(stem_object$dynamics$tcovar)
                        stem_object$dynamics$tcovar_adjmat       <- build_tcovar_adjmat(stem_object$dynamics$rates, stem_object$dynamics$tcovar_codes)
                }

                # generate or copy the initial states
                if(stem_object$dynamics$fixed_inits) {

                        if(stem_object$dynamics$n_strata == 1) {
                                param_inds     <- stem_object$dynamics$state_initializer$param_inds
                                initdist_codes <- stem_object$dynamics$state_initializer$codes
                        } else {
                                param_inds     <- unlist(lapply(stem_object$dynamics$state_initializer, function(x) x$param_inds))
                                initdist_codes <- unlist(lapply(stem_object$dynamics$state_initializer, function(x) x$codes))
                        }

                        # if all initial states are fixed, just copy the initial compartment counts
                        init_states    <- matrix(rep(stem_object$dynamics$initdist_params[param_inds], nsim), nrow = nsim, byrow = TRUE)[,order(initdist_codes), drop = FALSE]
                        colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                } else {
                        if(stem_object$dynamics$n_strata == 1) {

                                # simulate the initial compartment counts
                                init_states <- t(as.matrix(rmultinom(nsim, stem_object$dynamics$popsize, stem_object$dynamics$initdist_params)))
                                colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                        } else if(stem_object$dynamics$n_strata > 1) {

                                # generate the matrix of initial compartment counts
                                init_states <- matrix(0, nrow = nsim, ncol = length(stem_object$dynamics$comp_codes))
                                colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                                for(s in seq_len(stem_object$dynamics$n_strata)) {
                                        init_states[,stem_object$dynamics$state_initializer[[s]]$codes] <- as.matrix(t(rmultinom(nsim, stem_object$dynamics$strata_sizes[s], stem_object$dynamics$state_initializer[[s]]$init_states)))
                                }
                        }
                }

                # if there are artificial incidence compartments, copy the
                # incidence counts and add them to the initial state matrix
                if(!is.null(stem_object$dynamics$incidence_codes)) {
                        # init_incid <- init_states[, stem_object$dynamics$incidence_sources + 1, drop = FALSE]
                        init_incid <- matrix(0, nrow = nrow(init_states), ncol = length(stem_object$dynamics$incidence_codes))
                        colnames(init_incid) <- names(stem_object$dynamics$incidence_codes)
                        init_states <- cbind(init_states, init_incid)
                }

                # initialize the list of paths
                paths_full <- vector(mode = "list", length = nsim)

                # guess the initial dimensions. need an extra column for event times and another for event IDs.
                if(stem_object$dynamics$progressive & any(stem_object$dynamics$absorbing_states)) {
                        if(stem_object$dynamics$n_strata == 1) {
                                init_dims <- c(n_rows = sum(stem_object$dynamics$popsize * stem_object$dynamics$n_compartments),
                                               n_cols = ncol(stem_object$dynamics$flow_matrix) + 2)
                        } else {
                                init_dims <- c(n_rows = sum(stem_object$dynamics$strata_sizes * sapply(sapply(stem_object$dynamics$state_initializer, "[[", 4), length)),
                                               n_cols = ncol(stem_object$dynamics$flow_matrix) + 2)

                        }
                } else {
                        if(stem_object$dynamics$n_strata == 1) {
                                init_dims <- c(n_rows = sum(stem_object$dynamics$popsize * stem_object$dynamics$n_compartments) * 3,
                                               n_cols = ncol(stem_object$dynamics$flow_matrix) + 2)
                        } else if(stem_object$dynamics$n_strata > 1) {
                                init_dims <- c(n_rows = sum(stem_object$dynamics$strata_sizes * sapply(sapply(stem_object$dynamics$state_initializer, "[[", 4), length)) * 3,
                                               n_cols = ncol(stem_object$dynamics$flow_matrix) + 2)
                        }
                }

                # make the initial dimensions a little bigger (round up to nearest power of 2)
                p <- 1
                while(2^p < init_dims[1]) {
                        p <- p+1
                        if(2^p > init_dims[1]) init_dims[1] <- 2^p
                }

                # get the compartment names
                # path_colnames <- c("time", "event", c(names(stem_object$dynamics$comp_codes), names(stem_object$dynamics$incidence_codes)))

                # simulate the paths
                for(k in seq_len(nsim)) {
                        paths_full[[k]] <- simulate_gillespie(flow             = stem_object$dynamics$flow_matrix,
                                                              parameters       = stem_object$dynamics$parameters,
                                                              constants        = stem_object$dynamics$constants,
                                                              tcovar           = stem_object$dynamics$tcovar,
                                                              init_states      = init_states[k,],
                                                              rate_adjmat      = stem_object$dynamics$rate_adjmat,
                                                              tcovar_adjmat    = stem_object$dynamics$tcovar_adjmat,
                                                              tcovar_changemat = stem_object$dynamics$tcovar_changemat,
                                                              init_dims        = init_dims,
                                                              rate_ptr         = stem_object$dynamics$rate_ptrs[[1]])
                        # colnames(paths_full[[k]]) <- path_colnames
                }

                if(!is.null(census_times)) {
                        census_paths    <- vector(mode = "list", length = nsim)
                        # census_colnames <- c("time", c(names(stem_object$dynamics$comp_codes), names(stem_object$dynamics$incidence_codes)))

                        # add 2 to the codes b/c 'time' and 'event' are in the full path
                        census_codes      <- c(stem_object$dynamics$comp_codes, stem_object$dynamics$incidence_codes) + 2
                        get_incidence <- !is.null(stem_object$dynamics$incidence_codes)

                        if(get_incidence) {
                                incidence_codes       <- stem_object$dynamics$incidence_codes + 1
                                census_incidence_rows <- rep(list(seq_along(census_times) - 1), length(stem_object$dynamics$incidence_codes))
                        }


                        for(k in seq_len(nsim)) {
                                # get the census path
                                census_paths[[k]] <- build_census_path(path = paths_full[[k]], census_times = census_times, census_columns = census_codes)

                                # compute incidence if required. n.b. add 1 to the incidence codes b/c 'time' is in the census path
                                if(get_incidence) compute_incidence(censusmat = census_paths[[k]],
                                                                    col_inds  = incidence_codes,
                                                                    row_inds  = census_incidence_rows)

                                # assign column names
                                # colnames(census_paths[[k]]) <- census_colnames
                        }

                        if(paths_as_array) {
                                census_paths <- array(unlist(census_paths), dim = c(nrow(census_paths[[1]]), ncol(census_paths[[1]]), length(census_paths)))
                                # colnames(census_paths) <- census_colnames
                        }
                }

        } else if (method == "lna") {

                # set the restart times if they were not supplied
                if(identical(lna_restart, FALSE)) {
                        restart_times <- min(census_times)
                        lna_restart   <- FALSE
                } else {
                        if(is.logical(lna_restart)) {
                                restart_times <- census_times
                                lna_restart   <- TRUE
                        } else {
                                restart_times <- sort(unique(c(lna_restart, range(census_times))))
                                lna_restart   <- TRUE
                        }
                }

                # set the vectors of times when the LNA is evaluated, restarted, and censused
                lna_times       <- sort(unique(c(census_times, restart_times, stem_object$dynamics$.dynamics_args$tcovar[,1])))
                restart_inds    <- !is.na(match(lna_times, restart_times)) # indicator vector for when the LNA is to be restarted
                census_inds     <- !is.na(match(lna_times, census_times))

                # generate the matrix of parameters, constants, and time-varying covariates
                lna_pars <- matrix(0.0, nrow = length(lna_times), ncol = length(stem_object$dynamics$lna_param_codes))
                colnames(lna_pars) <- c(names(stem_object$dynamics$lna_param_codes))

                # insert parameters, constants, and time-varying covariates
                parameter_inds <- 1:length(stem_object$dynamics$param_codes)
                constant_inds  <- (length(parameter_inds)+1):(length(parameter_inds) + length(stem_object$dynamics$const_codes))
                tcovar_inds  <- (max(constant_inds)+1):ncol(lna_pars)
                lna_pars[, parameter_inds] <- matrix(stem_object$dynamics$parameters, nrow = nrow(lna_pars), ncol = length(parameter_inds), byrow = T)
                lna_pars[, constant_inds]  <- matrix(stem_object$dynamics$constants, nrow = nrow(lna_pars), ncol = length(constant_inds), byrow = T)

                if(!is.null(stem_object$dynamics$.dynamics_args$tcovar)) {
                        tcovar_rowinds          <- findInterval(lna_times, stem_object$dynamics$.dynamics_args$tcovar[,1])
                        lna_pars[, tcovar_inds] <- stem_object$dynamics$.dynamics_args$tcovar[tcovar_rowinds,-1]
                }

                # generate the initial values
                if(stem_object$dynamics$fixed_inits) {

                        if(stem_object$dynamics$n_strata == 1) {
                                param_inds     <- stem_object$dynamics$state_initializer$param_inds
                                initdist_codes <- stem_object$dynamics$state_initializer$codes
                        } else {
                                param_inds     <- unlist(lapply(stem_object$dynamics$state_initializer, function(x) x$param_inds))
                                initdist_codes <- unlist(lapply(stem_object$dynamics$state_initializer, function(x) x$codes))
                        }

                        # if all initial states are fixed, just copy the initial compartment counts
                        init_states <- matrix(rep(stem_object$dynamics$initdist_params[param_inds], nsim), nrow = nsim, byrow = TRUE)[,order(initdist_codes), drop = FALSE]
                        colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                } else {
                        if(stem_object$dynamics$n_strata == 1) {

                                # simulate the initial compartment counts
                                init_states <- t(as.matrix(rmultinom(nsim, stem_object$dynamics$popsize, stem_object$dynamics$initdist_params)))
                                colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                        } else if(stem_object$dynamics$n_strata > 1) {

                                # generate the matrix of initial compartment counts
                                init_states <- matrix(0, nrow = nsim, ncol = stem_object$dynamics$n_compartments)
                                colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                                for(s in seq_len(stem_object$dynamics$n_strata)) {
                                        init_states[,stem_object$dynamics$state_initializer[[s]]$codes] <- as.matrix(t(rmultinom(nsim, stem_object$dynamics$strata_sizes[s], stem_object$dynamics$state_initializer[[s]]$init_states)))
                                }
                        }
                }

                # if there are artificial incidence compartments, copy the
                # incidence counts and add them to the initial state matrix
                if(!is.null(stem_object$dynamics$incidence_codes)) {
                        init_incid <- init_states[, stem_object$dynamics$incidence_sources + 1, drop = FALSE]
                        colnames(init_incid) <- names(stem_object$dynamics$incidence_codes)
                        init_states <- cbind(init_states, init_incid)
                        lna_incidence_codes <- stem_object$dynamics$incidence_codes + 1 # add 1 b/c of the time column
                } else {
                        lna_incidence_codes <- 0
                }

                # if the LNA is on the log scale, log transform the initial states after adding a bit of noise
                if(stem_object$dynamics$lna_scale == "log") init_states <- log(init_states + 1e-6)

                census_paths <- simulate_lna(nsim             = nsim,
                                             lna_times        = lna_times,
                                             census_times     = census_times,
                                             census_inds      = census_inds,
                                             restart_inds     = restart_inds,
                                             lna_pars         = lna_pars,
                                             init_states      = init_states,
                                             incidence_codes  = lna_incidence_codes,
                                             log_scale        = stem_object$dynamics$lna_scale == "log",
                                             flow_matrix      = stem_object$dynamics$flow_matrix,
                                             lna_pointer      = stem_object$dynamics$lna_pointer$lna_ptr)

                colnames(census_paths) <- c("time", colnames(stem_object$dynamics$flow_matrix))
        }

        if(observations) {
                datasets         <- vector(mode = "list", length = nsim) # list for storing the datasets
                already_censused <- !is.null(census_times) && identical(stem_object$measurement_process$obstimes, census_times)
                measvar_names    <- colnames(stem_object$measurement_process$obsmat)

                # grab the time-varying covariate values at observation times
                tcovar_obstimes <- build_census_path(path           = stem_object$dynamics$tcovar,
                                                     census_times   = stem_object$measurement_process$obstimes,
                                                     census_columns = 1:(ncol(stem_object$dynamics$tcovar)-1))
                colnames(tcovar_obstimes) <- colnames(stem_object$dynamics$tcovar)

                if(already_censused) {

                        for(k in seq_len(nsim)) {
                                if(paths_as_array) {
                                        datasets[[k]] <- simulate_r_measure(censusmat = census_paths[,,k],
                                                                           measproc_indmat = stem_object$measurement_process$measproc_indmat,
                                                                           parameters = stem_object$dynamics$parameters,
                                                                           constants = stem_object$dynamics$constants,
                                                                           tcovar = tcovar_obstimes,
                                                                           r_measure_ptr = stem_object$measurement_process$meas_pointers$r_measure_ptr)
                                } else {
                                        datasets[[k]] <- simulate_r_measure(censusmat = census_paths[[k]],
                                                                           measproc_indmat = stem_object$measurement_process$measproc_indmat,
                                                                           parameters = stem_object$dynamics$parameters,
                                                                           constants = stem_object$dynamics$constants,
                                                                           tcovar = tcovar_obstimes,
                                                                           r_measure_ptr = stem_object$measurement_process$meas_pointers$r_measure_ptr)
                                }

                                colnames(datasets[[k]]) <- measvar_names
                        }


                } else if(!already_censused) {

                        censusmat       <- stem_object$measurement_process$censusmat
                        census_codes    <- c(stem_object$dynamics$comp_codes, stem_object$dynamics$incidence_codes) + 2

                        get_incidence <- !is.null(stem_object$dynamics$incidence_codes)

                        if(get_incidence) incidence_codes <- stem_object$dynamics$incidence_codes + 1

                        for(k in seq_len(nsim)) {
                                retrieve_census_path(censusmat, paths_full[[k]], stem_object$measurement_process$obstimes, census_columns = census_codes)
                                if(get_incidence) compute_incidence(censusmat = censusmat,
                                                                        col_inds  = incidence_codes,
                                                                        row_inds =  stem_object$measurement_process$obstime_inds)

                                datasets[[k]] <- simulate_r_measure(censusmat = censusmat,
                                                                   measproc_indmat = stem_object$measurement_process$measproc_indmat,
                                                                   parameters = stem_object$dynamics$parameters,
                                                                   constants = stem_object$dynamics$constants,
                                                                   tcovar = tcovar_obstimes,
                                                                   r_measure_ptr = stem_object$measurement_process$meas_pointers$r_measure_ptr)
                                colnames(datasets[[k]]) <- measvar_names
                        }
                }

                if(datasets_as_array) {
                        datasets <- array(unlist(datasets), dim = c(nrow(datasets[[1]]), ncol(datasets[[1]]), length(datasets)))
                        colnames(datasets) <- measvar_names
                }
        }

        stem_simulations <- list(paths = NULL, datasets = NULL, subject_paths = NULL)

        if(paths & is.null(census_times))  stem_simulations$paths <- paths_full
        if(paths & !is.null(census_times)) stem_simulations$paths <- census_paths
        if(observations)  stem_simulations$datasets      <- datasets
        if(subject_paths) stem_simulations$subject_paths <- subject_paths

        return(stem_simulations)
}