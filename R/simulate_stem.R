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
simulate_stem <-
        function(stem_object,
                 nsim = 1,
                 paths = FALSE,
                 observations = FALSE,
                 subject_paths = FALSE,
                 method = "gillespie",
                 tmax = NULL,
                 timestep = NULL,
                 census_times = NULL,
                 lna_path_type = "incidence",
                 messages = TRUE) {

                # ensure that the method is correctly specified
                if(!method %in% c("gillespie", "lna")) {
                        stop("The simulation method must either be 'gillespie' or 'lna'.")
                }

                if(method == "lna" && !lna_path_type %in% c("incidence", "prevalence", "both")) {
                        stop("The LNA path type must be specified as one of 'incidence', 'prevalence', or 'both'.")
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

                # t0 and tmax if not supplied
                t0                     <- stem_object$dynamics$t0
                if(is.null(tmax)) tmax <- stem_object$dynamics$tmax

                # make sure that t0 and tmax are in the census times
                if(!is.null(census_times)) census_times <- as.numeric(sort(unique(c(t0, census_times, tmax))))

                # build the time varying covariate matrix (includes, at a minimum, the endpoints of the simulation interval)
                # if timestep is null, there are no time-varying covariates
                if(method == "gillespie") {

                        # if any of t0, tmax, or a timestep was supplied,
                        # check if they differ from the parameters supplied in the stem_object$dynamics.
                        # if they differ, reconstruct the tcovar matrix and associated objects
                        rebuild_tcovar <- (t0 != stem_object$dynamics$t0) ||
                                (tmax != stem_object$dynamics$tcovar[nrow(stem_object$dynamics$tcovar), 1]) ||
                                !is.null(timestep)

                        if(rebuild_tcovar) {

                                # replace t0 and tmax in the time-varying covariate matrix
                                if(is.null(timestep) && stem_object$dynamics$timevarying) {
                                        timestep <- (tmax - t0)/50

                                } else if(is.null(timestep) && !stem_object$dynamics$timevarying){
                                        timestep <- tmax - t0
                                }

                                # rebuild the time-varying covariate matrix so that it contains the census intervals
                                stem_object$dynamics$tcovar <- build_tcovar_matrix(tcovar = stem_object$dynamics$.dynamics_args$tcovar,
                                                                                   timestep = timestep, t0 = t0, tmax = tmax)
                                stem_object$dynamics$tcovar_codes        <- seq_len(ncol(stem_object$dynamics$tcovar) - 1)
                                names(stem_object$dynamics$tcovar_codes) <- colnames(stem_object$dynamics$tcovar)[2:ncol(stem_object$dynamics$tcovar)]
                                stem_object$dynamics$n_tcovar            <- ncol(stem_object$dynamics$tcovar) - 1
                                stem_object$dynamics$tcovar_changemat    <- build_tcovar_changemat(stem_object$dynamics$tcovar)
                                stem_object$dynamics$tcovar_adjmat       <- build_tcovar_adjmat(stem_object$dynamics$rates, stem_object$dynamics$tcovar_codes)
                        }

                        # generate or copy the initial states
                        if(stem_object$dynamics$fixed_inits) {

                                # if all initial states are fixed, just copy the initial compartment counts
                                init_states <- matrix(
                                        rep(stem_object$dynamics$initdist_params, nsim),
                                        nrow = nsim,
                                        byrow = TRUE
                                )
                                colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                        } else {

                                if(stem_object$dynamics$n_strata == 1) {

                                        # simulate the initial compartment counts
                                        init_states <- t(as.matrix(rmultinom(nsim, stem_object$dynamics$popsize,stem_object$dynamics$initdist_priors)))
                                        colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                                } else if(stem_object$dynamics$n_strata > 1) {

                                        # generate the matrix of initial compartment counts
                                        init_states <- matrix(0, nrow = nsim, ncol = length(stem_object$dynamics$comp_codes))
                                        colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                                        for(s in seq_len(stem_object$dynamics$n_strata)) {
                                                init_states[,stem_object$dynamics$state_initializer[[s]]$codes] <- as.matrix(t(rmultinom(nsim, stem_object$dynamics$strata_sizes[s], stem_object$dynamics$state_initializer[[s]]$prior)))
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
                        path_colnames <-
                                c("time", "event", c(
                                        names(stem_object$dynamics$comp_codes),
                                        names(stem_object$dynamics$incidence_codes)
                                ))

                        if(is.null(census_times)) {
                                # initialize the list of paths
                                paths_full <- vector(mode = "list", length = nsim)

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
                                        colnames(paths_full[[k]]) <- path_colnames
                                }

                        } else {
                                census_paths    <- vector(mode = "list", length = nsim)
                                census_colnames <- c("time", c(names(stem_object$dynamics$comp_codes),
                                                               names(stem_object$dynamics$incidence_codes)))

                                # add 2 to the codes b/c 'time' and 'event' are in the full path
                                census_codes  <- c(stem_object$dynamics$comp_codes, stem_object$dynamics$incidence_codes) + 2
                                get_incidence <- !is.null(stem_object$dynamics$incidence_codes)

                                if(get_incidence) {
                                        incidence_codes       <- stem_object$dynamics$incidence_codes + 1
                                        census_incidence_rows <- rep(list(seq_along(census_times) - 1),
                                                                     length(stem_object$dynamics$incidence_codes))
                                }

                                for(k in seq_len(nsim)) {

                                        path_full <- simulate_gillespie(flow             = stem_object$dynamics$flow_matrix,
                                                                        parameters       = stem_object$dynamics$parameters,
                                                                        constants        = stem_object$dynamics$constants,
                                                                        tcovar           = stem_object$dynamics$tcovar,
                                                                        init_states      = init_states[k,],
                                                                        rate_adjmat      = stem_object$dynamics$rate_adjmat,
                                                                        tcovar_adjmat    = stem_object$dynamics$tcovar_adjmat,
                                                                        tcovar_changemat = stem_object$dynamics$tcovar_changemat,
                                                                        init_dims        = init_dims,
                                                                        rate_ptr         = stem_object$dynamics$rate_ptrs[[1]])

                                        # get the census path
                                        census_paths[[k]] <- build_census_path(path = path_full,
                                                                               census_times = census_times,
                                                                               census_columns = census_codes)

                                        # compute incidence if required. n.b. add 1 to the incidence codes b/c 'time' is in the census path
                                        if(get_incidence) compute_incidence(censusmat = census_paths[[k]],
                                                                            col_inds  = incidence_codes,
                                                                            row_inds  = census_incidence_rows)

                                        # assign column names
                                        colnames(census_paths[[k]]) <- census_colnames
                                }
                        }

                } else if (method == "lna") {

                        # set the vectors of times when the LNA is evaluated and censused
                        lna_times <- sort(unique(c(t0, census_times, stem_object$dynamics$.dynamics_args$tcovar[,1], tmax)))

                        # generate the matrix of parameters, constants, and time-varying covariates
                        lna_pars  <- matrix(0.0,
                                            nrow = length(lna_times),
                                            ncol = length(stem_object$dynamics$lna_rates$lna_param_codes),
                                            dimnames = list(NULL, names(stem_object$dynamics$lna_rates$lna_param_codes)))

                        # insert parameters, constants, and time-varying covariates
                        parameter_inds <- 1:length(stem_object$dynamics$param_codes)
                        constant_inds  <- (length(parameter_inds)+1):(length(parameter_inds) + length(stem_object$dynamics$const_codes))
                        tcovar_inds    <- (max(constant_inds)+1):ncol(lna_pars)
                        lna_pars[, parameter_inds] <- matrix(stem_object$dynamics$parameters,
                                                             nrow = nrow(lna_pars),
                                                             ncol = length(parameter_inds), byrow = T)
                        lna_pars[, constant_inds]  <- matrix(stem_object$dynamics$constants,
                                                             nrow = nrow(lna_pars),
                                                             ncol = length(constant_inds), byrow = T)

                        if(!is.null(stem_object$dynamics$.dynamics_args$tcovar)) {
                                tcovar_rowinds          <- findInterval(lna_times, stem_object$dynamics$.dynamics_args$tcovar[,1])
                                lna_pars[, tcovar_inds] <- stem_object$dynamics$.dynamics_args$tcovar[tcovar_rowinds,-1]
                        }

                        # generate some auxilliary objects
                        param_update_inds <- is.na(match(lna_times, census_times))

                        # simulate the LNA paths
                        census_paths <- vector(mode = "list", length = nsim)
                        lna_paths    <- vector(mode = "list", length = nsim)

                        # retrieve the initial state
                        init_state <- lna_pars[1, stem_object$dynamics$lna_initdist_inds + 1]

                        # if the initial state is not fixed, sample the collection of initial states
                        if(!stem_object$dynamics$fixed_inits) {
                                if(stem_object$dynamics$n_strata == 1) {

                                        # simulate the initial compartment counts
                                        init_states <- stem_object$dynamics$popsize * extraDistr::rdirichlet(nsim, stem_object$dynamics$initdist_priors)
                                        colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                                } else if(stem_object$dynamics$n_strata > 1) {

                                        # generate the matrix of initial compartment counts
                                        init_states <- matrix(0.0, nrow = nsim, ncol = length(stem_object$dynamics$comp_codes))
                                        colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                                        for(s in seq_len(stem_object$dynamics$n_strata)) {
                                                init_states[,stem_object$dynamics$state_initializer[[s]]$codes] <-
                                                        stem_object$dynamics$strata_sizes[s] * extraDistr::rdirichlet(nsim, stem_object$dynamics$state_initializer[[s]]$prior)
                                        }
                                }
                        }

                        for(k in seq_len(nsim)) {

                                if(!stem_object$dynamics$fixed_inits) {
                                        pars2lnapars(lna_pars, c(parameters, init_states[k,]))
                                }

                                census_paths[[k]] <- propose_lna(lna_times         = lna_times,
                                                                 lna_pars          = lna_pars,
                                                                 init_start        = stem_object$dynamics$lna_initdist_inds[1],
                                                                 param_update_inds = param_update_inds,
                                                                 stoich_matrix     = stem_object$dynamics$stoich_matrix_lna,
                                                                 net_effect        = stem_object$dynamics$net_effect,
                                                                 lna_pointer       = stem_object$dynamics$lna_pointers$lna_ptr,
                                                                 set_pars_pointer  = stem_object$dynamics$lna_pointers$set_lna_params_ptr)$lna_path

                                lna_paths[[k]] <- lna_incid2prev(census_paths[[k]],
                                                              stem_object$dynamics$flow_matrix_lna,
                                                              init_states[k,])

                                colnames(census_paths[[k]]) <- c("time", rownames(stem_object$dynamics$flow_matrix_lna))
                                colnames(lna_paths[[k]]) <- c("time",colnames(stem_object$dynamics$flow_matrix_lna))
                        }
                }

                if(observations) {

                        datasets         <- vector(mode = "list", length = nsim) # list for storing the datasets
                        measvar_names    <- colnames(stem_object$measurement_process$obsmat)

                        # grab the time-varying covariate values at observation times
                        tcovar_obstimes <- build_census_path(path           = stem_object$dynamics$tcovar,
                                                             census_times   = stem_object$measurement_process$obstimes,
                                                             census_columns = 1:(ncol(stem_object$dynamics$tcovar)-1))
                        colnames(tcovar_obstimes) <- colnames(stem_object$dynamics$tcovar)

                        # if incidence, the incidence codes are not null
                        do_incidence <- !is.null(stem_object$dynamics$incidence_codes) & method != "lna"

                        # census if computing incidence or if computing prevalence and the obstimes != census_times
                        do_census <- do_incidence ||
                                     (!is.null(census_times) &&
                                      !identical(as.numeric(stem_object$measurement_process$obstimes), as.numeric(census_times)))

                        # get the indices in the censused matrices for the observation times
                        if(do_census) cens_inds <- findInterval(stem_object$measurement_process$obstimes, census_times)

                        if(method == "gillespie") {

                                # should incidence be computed?
                                if(do_incidence) incidence_codes <- stem_object$dynamics$incidence_codes + 1

                                # column codes in the path matrix for compartments to be censused
                                census_codes    <- c(stem_object$dynamics$comp_codes,
                                                     stem_object$dynamics$incidence_codes) + 2
                                censusmat       <- stem_object$measurement_process$censusmat

                                for(k in seq_len(nsim)) {

                                        if(is.null(census_times)) {
                                                # get the state at the observation times
                                                cens_inds <- findInterval(stem_object$measurement_process$obstimes,
                                                                          paths_full[[k]][,1])

                                                # fill out the census matrix
                                                censusmat[,-1] <- paths_full[[k]][cens_inds, census_codes + 1]

                                                # compute the incidence if appropriate
                                                if(get_incidence) compute_incidence(censusmat = censusmat,
                                                                                    col_inds  = incidence_codes,
                                                                                    row_inds  = stem_object$measurement_process$obstime_inds)
                                        } else {
                                                # get the state at observation times
                                                if(do_census) {
                                                        censusmat[,-1] <- census_paths[[k]][cens_inds, -1]
                                                } else {
                                                        censusmat <- census_paths[[k]]
                                                }
                                        }

                                        # simulate the data
                                        datasets[[k]] <- simulate_r_measure(censusmat = censusmat,
                                                                            measproc_indmat = stem_object$measurement_process$measproc_indmat,
                                                                            parameters = stem_object$dynamics$parameters,
                                                                            constants = stem_object$dynamics$constants,
                                                                            tcovar = tcovar_obstimes,
                                                                            r_measure_ptr = stem_object$measurement_process$meas_pointers$r_measure_ptr)
                                        colnames(datasets[[k]]) <- measvar_names
                                }

                        } else if(method == "lna") {

                                # get the objects for simulating from the measurement process
                                measproc_indmat  <- stem_object$measurement_process$measproc_indmat
                                parameters       <- stem_object$dynamics$parameters
                                constants        <- stem_object$dynamics$parameters
                                tcovar           <- tcovar_obstimes
                                r_measure_ptr    <- stem_object$measurement_process$meas_pointers$r_measure_ptr
                                cens_inds        <- findInterval(stem_object$measurement_process$obstimes, census_times) - 1
                                do_prevalence    <- stem_object$measurement_process$lna_prevalence
                                do_incidence     <- stem_object$measurement_process$lna_incidence
                                incidence_codes  <- stem_object$measurement_process$incidence_codes_lna
                                obstime_inds     <- stem_object$measurement_process$obstime_inds
                                pathmat          <- stem_object$measurement_process$censusmat
                                flow_matrix_lna  <- stem_object$dynamics$flow_matrix_lna
                                cens_incid_codes <- incidence_codes + ncol(pathmat) - (1 + length(incidence_codes))

                                for(k in seq_len(nsim)) {

                                        # fill out the census matrix
                                        census_lna(path                = census_paths[[k]],
                                                   census_path         = pathmat,
                                                   census_inds         = cens_inds,
                                                   flow_matrix_lna     = flow_matrix_lna,
                                                   do_prevalence       = do_prevalence,
                                                   init_state          = init_state,
                                                   incidence_codes_lna = incidence_codes)
#
                                        # simulate the dataset
                                        datasets[[k]] <- simulate_r_measure(pathmat,
                                                                            measproc_indmat,
                                                                            parameters,
                                                                            constants,
                                                                            tcovar_obstimes,
                                                                            r_measure_ptr)
                                        colnames(datasets[[k]]) <- measvar_names
                                }
                        }
                }

                stem_simulations <- list(paths = NULL, datasets = NULL, subject_paths = NULL)

                if(paths & is.null(census_times)) {
                        stem_simulations$paths <- paths_full
                } else if(paths & !is.null(census_times)) {
                        stem_simulations$paths <- census_paths
                        if(method == "lna") stem_simulations$natural_paths <- lna_paths
                } else {
                        stem_simulations$paths <- NULL
                }
                if(observations)  stem_simulations$datasets      <- datasets
                if(subject_paths) stem_simulations$subject_paths <- subject_paths

                return(stem_simulations)
        }