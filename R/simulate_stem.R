#' Simulations from a stochastic epidemic model.
#'
#' @param stem Stochastic epidemic model object. The dynamics slot must be
#'   supplied at a minimum. If simulations from the measurement process are
#'   desired, the measurement process slot must also be supplied.
#' @param nsim number of realizations to simulate
#' @param simulation_parameters optional list of vectors of simulation
#'   parameters. If NULL, paths are simulated using the parameters specified in
#'   the stem object.
#' @param paths Should population-level paths at census times be returned?
#' @param full_paths Should complete population level paths be returned (only
#'   for Gillespie), defaults to FALSE
#' @param observations Should simulated observations be returned? Requires that
#'   a measurement process be defined in the stem object.
#' @param subject_paths Should population-level paths be mapped to subject-level
#'   paths and returned (e.g. for the purpose of initializing a subject-level
#'   collection of disease histories)? Only available for exact simulation via
#'   the Gillespie direct method.
#' @param method either "gillespie" if simulating via Gillespie's direct method,
#'   "lna" if simulating paths via the linear noise approximation, or "ode" if
#'   simulating paths of the deterministic limit of the underlying Markov jump
#'   process.
#' @param tmax the time at which simulation of the system is terminated. If not
#'   supplied, defaults to the last observation time if not supplied.
#' @param census_times vector of times at which compartment counts should be
#'   recorded.
#' @param max_attempts maximum number of times to attempt simulating an LNA path
#'   before aborting due to the moments of the path being degenerate at some
#'   point
#' @param lna_method defaults to "exact". If "approx", an initial path is
#'   proposed by resampling perturbations that lead to negative increments or
#'   volumes. The initial path is then updated via elliptical slice sampling
#' @param ess_warmup number of elliptical slice sampling updates before the
#'   lna sample is saved
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
                 simulation_parameters = NULL,
                 paths = FALSE,
                 full_paths = FALSE,
                 observations = FALSE,
                 subject_paths = FALSE,
                 method = "gillespie",
                 tmax = NULL, 
                 census_times = NULL,
                 max_attempts = 500,
                 lna_method = "exact",
                 ess_updates = 10,
                 ess_warmup = 100,
                 messages = TRUE) {

                # ensure that the method is correctly specified
                if(!method %in% c("gillespie", "lna", "ode")) {
                        stop("The simulation method must either be 'gillespie', 'lna', or 'ode'.")
                }

                # if lna, subject paths are not available
                if(subject_paths && (method == "lna" | method == "ode")) {
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

                # check that if simulation parameters are supplied, they are supplied in a list of named vectors of
                # length equal to the number of simulations. Also check that the parameters are supplied in the same
                # order as the parameters in the stem object.
                if(!is.null(simulation_parameters)) {
                        if(!is.list(simulation_parameters)) {
                                stop("Simulation parameters must be supplied as a list of vectors.")
                        }

                        if(length(simulation_parameters) != nsim) {
                                stop("If simulation parameters are supplied, there must be an equal number of vectors in the list as requested paths to be simulated.")
                        }

                        if(!all(sapply(simulation_parameters, function(x) all(names(x) == names(stem_object$dynamics$parameters))))) {
                                stop("The simulation parameters list must consist of named vectors with elements given in the same order as the parameters in the stem object.")
                        }
                }

                if(paths == observations) {
                        paths <- TRUE
                        observations <- TRUE
                }

                # t0 and tmax if not supplied
                t0                     <- stem_object$dynamics$t0
                if(is.null(tmax)) tmax <- stem_object$dynamics$tmax

                if(is.null(stem_object$dynamics$timestep)) {
                        timestep <- 1
                } else {
                        timestep <- stem_object$dynamics$timestep
                }

                # make sure that there exists a vector of census times
                if(is.null(census_times)) {
                        if(is.null(stem_object$measurement_process$obstimes)) {
                                census_times <- as.numeric(unique(c(t0, seq(t0,tmax,by=timestep), tmax)))

                        } else {
                                census_times <- as.numeric(unique(c(t0, stem_object$measurement_process$obstimes, tmax)))
                        }

                } else {
                        census_times <- as.numeric(sort(unique(c(t0, census_times, tmax))))
                }

                # build the time varying covariate matrix (includes, at a minimum, the endpoints of the simulation interval)
                # if timestep is null, there are no time-varying covariates
                if(method == "gillespie") {

                        # if any of t0, tmax, or a timestep was supplied,
                        # check if they differ from the parameters supplied in the stem_object$dynamics.
                        # if they differ, reconstruct the tcovar matrix and associated objects
                        rebuild_tcovar <- (t0 != stem_object$dynamics$t0) ||
                                (tmax != stem_object$dynamics$tcovar[nrow(stem_object$dynamics$tcovar), 1]) ||
                                timestep != stem_object$dynamics$timestep ||
                                !identical(stem_object$dynamics$tcovar[,1], census_times)

                        if(rebuild_tcovar) {

                                # rebuild the time-varying covariate matrix so that it contains the census intervals
                                stem_object$dynamics$tcovar <- build_tcovar_matrix(tcovar       = stem_object$dynamics$dynamics_args$tcovar, 
                                                                                   tparam       = stem_object$dynamics$tparam,
                                                                                   forcings     = stem_object$dynamics$forcings,
                                                                                   parameters   = stem_object$dynamics$parameters,
                                                                                   timestep     = timestep, 
                                                                                   census_times = census_times,
                                                                                   t0           = t0, 
                                                                                   tmax         = tmax, 
                                                                                   messages     = messages)
                                
                                stem_object$dynamics$tcovar_codes        <- seq_len(ncol(stem_object$dynamics$tcovar) - 1)
                                names(stem_object$dynamics$tcovar_codes) <- colnames(stem_object$dynamics$tcovar)[2:ncol(stem_object$dynamics$tcovar)]
                                stem_object$dynamics$n_tcovar            <- ncol(stem_object$dynamics$tcovar) - 1
                                stem_object$dynamics$tcovar_changemat    <- build_tcovar_changemat(tcovar = stem_object$dynamics$tcovar,
                                                                                                   tparam = stem_object$dynamics$tparam,
                                                                                                   forcings = stem_object$dynamics$forcings)
                                stem_object$dynamics$tcovar_adjmat       <- build_tcovar_adjmat(rates        = stem_object$dynamics$rates, 
                                                                                                tcovar_codes = stem_object$dynamics$tcovar_codes,
                                                                                                forcings     = stem_object$dynamics$forcings)

                                # zero out forcings if necessary
                                if(!is.null(stem_object$dynamics$dynamics_args$forcings)) {

                                        # get the forcing indices (supplied in the original tcovar matrix)
                                        forcing_inds <- vector("logical", nrow(stem_object$dynamics$tcovar))
                                        for(f in seq_along(stem_object$dynamics$forcings)) {
                                              forcing_inds <- forcing_inds | stem_object$dynamics$tcovar[,stem_object$dynamics$forcings[[f]]$tcovar_name] != 0
                                        }
                                        
                                } else {
                                        forcing_inds   <- rep(FALSE, nrow(stem_object$dynamics$tcovar))
                                }
                                
                        } else {
                                # Get the forcing indices if there are forcings
                                if(!is.null(stem_object$dynamics$dynamics_args$forcings)) {

                                        # get the forcing indices (supplied in the original tcovar matrix)
                                        forcing_inds <- vector("logical", nrow(stem_object$dynamics$tcovar))
                                        for(f in seq_along(stem_object$dynamics$forcings)) {
                                              forcing_inds <- forcing_inds | stem_object$dynamics$tcovar[,stem_object$dynamics$forcings[[f]]$tcovar_name] != 0
                                        }
                                        
                                } else {
                                        forcing_inds   <- rep(FALSE, nrow(stem_object$dynamics$tcovar))
                                }
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

                                                if(!stem_object$dynamics$state_initializer[[s]]$fixed) {
                                                        init_states[, stem_object$dynamics$state_initializer[[s]]$codes] <-
                                                                as.matrix(t(rmultinom(nsim,
                                                                                      stem_object$dynamics$strata_sizes[s],
                                                                                      stem_object$dynamics$state_initializer[[s]]$prior)))

                                                } else {
                                                        init_states[, stem_object$dynamics$state_initializer[[s]]$codes] <-
                                                                matrix(stem_object$dynamics$state_initializer[[s]]$init_states,
                                                                       nrow = nsim,
                                                                       ncol = length(stem_object$dynamics$state_initializer[[s]]$init_states),
                                                                       byrow = T)
                                                }
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

                        # generate forcings if they are specified
                        forcing_matrix <- matrix(0.0,
                                                 nrow = nrow(stem_object$dynamics$tcovar),
                                                 ncol = length(path_colnames)-2,
                                                 dimnames = list(NULL, path_colnames[-c(1:2)]))

                        if(!is.null(stem_object$dynamics$dynamics_args$forcings)) {

                                for(l in seq_along(stem_object$dynamics$dynamics_args$forcings)) {

                                        # insert the flow into the forcing matrix
                                        forcing_matrix[forcing_inds, stem_object$dynamics$forcings[[l]]$from] <-
                                                forcing_matrix[forcing_inds, stem_object$dynamics$forcings[[l]]$from] -
                                                stem_object$dynamics$tcovar[forcing_inds, stem_object$dynamics$forcings[[l]]$tcovar_name]

                                        forcing_matrix[forcing_inds, stem_object$dynamics$forcings[[l]]$to] <-
                                              forcing_matrix[forcing_inds, stem_object$dynamics$forcings[[l]]$to] +
                                              stem_object$dynamics$tcovar[forcing_inds, stem_object$dynamics$forcings[[l]]$tcovar_name]
                                }
                        }
                        
                        # tparam indices and initial values
                        if(!is.null(stem_object$dynamics$tparam)) {
                              
                              # generate indices for time-varying parameters
                              for(s in seq_along(stem_object$dynamics$tparam)) {
                                    stem_object$dynamics$tparam[[s]]$col_ind  <- 
                                          stem_object$dynamics$tcovar_codes[stem_object$dynamics$tparam[[s]]$tparam_name]
                                    stem_object$dynamics$tparam[[s]]$tpar_inds <- 
                                          findInterval(stem_object$dynamics$tcovar[,1], tparam[[s]]$times, left.open = F) - 1
                              }
                              
                              # list for saving the time varying parameters for reuse in simulating a dataset in necessary
                              tparam_draws <- vector("list", length = nsim)
                              tparam_times <- sort(unique(unlist(lapply(tparam, function(x) x$times))))
                              tparam_times <- tparam_times[tparam_times >= t0 & tparam_times <= tmax]
                        } else {
                              tparam_draws <- NULL
                              tparam_times <- NULL
                        }

                        # initialize the list of paths
                        if(full_paths) paths_full <- vector(mode = "list", length = nsim)

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

                        # vector of model parameters
                        sim_pars <- stem_object$dynamics$parameters

                        for(k in seq_len(nsim)) {

                                attempt <- 0
                                path_full <- NULL
                                if(!is.null(simulation_parameters)) {
                                      sim_pars <- as.numeric(simulation_parameters[[k]])
                                } 
                                
                                # draw new time-varying parameters if necessary
                                if(!is.null(stem_object$dynamics$tparam)) {
                                      for(s in seq_along(stem_object$dynamics$tparam)) {
                                            
                                            #draw new values
                                            draw_normals(stem_object$dynamics$tparam[[s]]$values)
                                            
                                            # insert the new values into the tcovar matrix
                                            insert_tparam(tcovar = stem_object$dynamics$tcovar, 
                                                          values = stem_object$dynamics$tparam[[s]]$draws2par(
                                                                sim_pars, 
                                                                stem_object$dynamics$tparam[[s]]$values),
                                                          col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                                                          tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds)
                                            
                                      }
                                      
                                      tparam_draws[[k]] <- lapply(stem_object$dynamics$tparam, function(x) x$values)
                                }

                                while(is.null(path_full) & attempt < max_attempts) {
                                        try({
                                                path_full <- simulate_gillespie(flow             = stem_object$dynamics$flow_matrix,
                                                                                parameters       = sim_pars,
                                                                                constants        = stem_object$dynamics$constants,
                                                                                tcovar           = stem_object$dynamics$tcovar,
                                                                                init_states      = init_states[k,],
                                                                                rate_adjmat      = stem_object$dynamics$rate_adjmat,
                                                                                tcovar_adjmat    = stem_object$dynamics$tcovar_adjmat,
                                                                                tcovar_changemat = stem_object$dynamics$tcovar_changemat,
                                                                                init_dims        = init_dims,
                                                                                forcing_inds     = forcing_inds,
                                                                                forcing_matrix   = forcing_matrix,
                                                                                rate_ptr         = stem_object$dynamics$rate_ptrs[[1]])
                                        }, silent = TRUE)
                                        attempt <- attempt + 1
                                }

                                # get the census path
                                if(!is.null(path_full)) {
                                      census_paths[[k]] <- build_census_path(path = path_full,
                                                                             census_times = census_times,
                                                                             census_columns = census_codes)
                                      
                                      # compute incidence if required. n.b. add 1 to the incidence codes b/c 'time' is in the census path
                                      if(get_incidence) compute_incidence(censusmat = census_paths[[k]],
                                                                          col_inds  = incidence_codes,
                                                                          row_inds  = census_incidence_rows)
                                      
                                      # assign column names
                                      colnames(census_paths[[k]]) <- census_colnames
                                      
                                      # save the full path if desired
                                      if(full_paths) paths_full[[k]] <- path_full
                                }
                        }
                        
                        failed_runs  <- which(sapply(census_paths, is.null))
                        if(length(failed_runs) != 0) {
                              census_paths <- census_paths[-failed_runs]
                              if(full_paths) paths_full <- paths_full[-failed_runs]
                        }

                } else if (method == "lna") {

                        # pull out logical vector for which rates are of order > 1
                        higher_order_rates <- sapply(stem_object$dynamics$rates, function(x) x$higher_order)

                        # set the vectors of times when the LNA is evaluated and censused
                        lna_times <- sort(unique(
                                        c(t0,
                                          census_times,
                                          seq(t0, tmax, by = stem_object$dynamics$timestep),
                                          stem_object$dynamics$tcovar[, 1],
                                          tmax)))

                        # generate the matrix of parameters, constants, and time-varying covariates
                        lna_pars  <- matrix(0.0,
                                            nrow = length(lna_times),
                                            ncol = length(stem_object$dynamics$lna_rates$lna_param_codes),
                                            dimnames = list(NULL, names(stem_object$dynamics$lna_rates$lna_param_codes)))

                        # insert parameters, constants, and time-varying covariates
                        parameter_inds <- seq_along(stem_object$dynamics$param_codes)
                        constant_inds  <- (length(parameter_inds)+1):(length(parameter_inds) + length(stem_object$dynamics$const_codes))
                        tcovar_inds    <- (max(constant_inds)+1):ncol(lna_pars)
                        lna_pars[, parameter_inds] <- matrix(stem_object$dynamics$parameters,
                                                             nrow = nrow(lna_pars),
                                                             ncol = length(parameter_inds), byrow = T)
                        lna_pars[, constant_inds]  <- matrix(stem_object$dynamics$constants,
                                                             nrow = nrow(lna_pars),
                                                             ncol = length(constant_inds), byrow = T)

                        # Generate forcing indices and update the LNA parameter matrix with forcings
                        forcing_inds <- rep(FALSE, length(lna_times))

                        if(!is.null(stem_object$dynamics$dynamics_args$tcovar)) {
                                tcovar_rowinds <- match(round(lna_times, digits = 8),
                                                        round(stem_object$dynamics$tcovar[,1], digits = 8))
                                lna_pars[tcovar_rowinds, tcovar_inds] <- stem_object$dynamics$tcovar[tcovar_rowinds,-1]

                                # zero out forcings if necessary
                                if(!is.null(stem_object$dynamics$dynamics_args$forcings)) {

                                        # get the forcing indices (supplied in the original tcovar matrix)
                                        forcing_inds <- as.logical(match(round(lna_times, digits = 8),
                                                                         round(stem_object$dynamics$dynamics_args$tcovar[,1], digits = 8),
                                                                         nomatch = FALSE))
                                        zero_inds    <- !forcing_inds

                                        # zero out the tcovar elements corresponding to times with no forcings
                                        for(l in seq_along(stem_object$dynamics$dynamics_args$forcings)) {
                                                lna_pars[zero_inds, stem_object$dynamics$dynamics_args$forcings[[l]]$tcovar_name] = 0
                                        }
                                }
                        }
                        
                        # tparam indices and initial values
                        if(!is.null(stem_object$dynamics$tparam)) {
                              
                              # generate indices for time-varying parameters
                              for(s in seq_along(stem_object$dynamics$tparam)) {
                                    stem_object$dynamics$tparam[[s]]$col_ind  <- 
                                          stem_object$dynamics$lna_rates$lna_param_codes[stem_object$dynamics$tparam[[s]]$tparam_name]
                                    stem_object$dynamics$tparam[[s]]$tpar_inds <- 
                                          findInterval(lna_times, tparam[[s]]$times, left.open = F) - 1
                              }
                              
                              # list for saving the time varying parameters for reuse in simulating a dataset in necessary
                              tparam_draws <- vector("list", length = nsim)
                              tparam_times <- sort(unique(unlist(lapply(tparam, function(x) x$times))))
                              tparam_times <- tparam_times[tparam_times >= t0 & tparam_times <= tmax]
                        } else {
                              tparam_draws <- NULL
                              tparam_times <- NULL
                        }

                        # generate some auxilliary objects
                        param_update_inds <- lna_times %in% 
                              sort(unique(c(t0, tmax, stem_object$dynamics$tcovar[,1], tparam_times))) 
                        census_interval_inds <- findInterval(lna_times, census_times, left.open = T)

                        # simulate the LNA paths
                        census_paths <- vector(mode = "list", length = nsim)
                        lna_draws    <- vector(mode = "list", length = nsim)
                        lna_paths    <- vector(mode = "list", length = nsim)
                        # drift_vecs   <- vector(mode = "list", length = nsim)
                        # diff_mats    <- vector(mode = "list", length = nsim)

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

                                                if(!stem_object$dynamics$state_initializer[[s]]$fixed) {

                                                        init_states[,stem_object$dynamics$state_initializer[[s]]$codes] <-
                                                                stem_object$dynamics$strata_sizes[s] * extraDistr::rdirichlet(nsim, stem_object$dynamics$state_initializer[[s]]$prior)

                                                } else {
                                                        init_states[,stem_object$dynamics$state_initializer[[s]]$codes] <-
                                                                matrix(stem_object$dynamics$state_initializer[[s]]$init_states,
                                                                       nrow = nsim,
                                                                       ncol = length(stem_object$dynamics$state_initializer[[s]]$init_states),
                                                                       byrow = TRUE)
                                                }
                                        }
                                }
                        } else {
                                if(stem_object$dynamics$n_strata == 1) {
                                        init_states <- matrix(rep(stem_object$dynamics$state_initializer$init_states, nsim), nrow = nsim, byrow = T)

                                } else {
                                        init_states <- matrix(0.0, nrow = nsim, ncol = length(stem_object$dynamics$comp_codes))
                                        colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                                        for(s in seq_len(stem_object$dynamics$n_strata)) {
                                                init_states[,stem_object$dynamics$state_initializer[[s]]$codes] <-
                                                        matrix(stem_object$dynamics$state_initializer[[s]]$init_states,
                                                               nrow = nsim,
                                                               ncol = length(stem_object$dynamics$state_initializer[[s]]$init_states),
                                                               byrow = TRUE)
                                        }
                                }
                        }

                        # generate forcings if they are specified
                        forcing_matrix <- matrix(0.0,
                                                 nrow = length(lna_times),
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

                        sim_pars  <- stem_object$dynamics$parameters
                        init_vols <- init_states[1,]
                        prev_inds <- match(round(census_times, digits = 8), round(lna_times, digits = 8))
                        
                        for(k in seq_len(nsim)) {
                              
                              if(!is.null(simulation_parameters)) {
                                    sim_pars <- as.numeric(simulation_parameters[[k]])
                              }
                              
                              if(!stem_object$dynamics$fixed_inits) {
                                    init_vols <- init_states[k,]
                              }
                              
                              # set the parameters and initial volumes
                              pars2lnapars(lna_pars, c(sim_pars, init_vols))
                              
                              # draw values for the time-varying parameters
                              if(!is.null(stem_object$dynamics$tparam)) {
                                    for(s in seq_along(stem_object$dynamics$tparam)) {
                                          
                                          #draw new values
                                          draw_normals(stem_object$dynamics$tparam[[s]]$values)
                                          
                                          # insert the new values into the tcovar matrix
                                          insert_tparam(tcovar    = lna_pars, 
                                                        values    = stem_object$dynamics$tparam[[s]]$draws2par(stem_object$dynamics$parameters,
                                                                                                               stem_object$dynamics$tparam[[s]]$values),
                                                        col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                                                        tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds)
                                    }
                                    
                                    tparam_draws[[k]] <- lapply(stem_object$dynamics$tparam, function(x) x$values)
                              }
                              
                              attempt <- 0
                              path    <- NULL
                              while(is.null(path) && (attempt < max_attempts)) {
                                    
                                    if(lna_method == "exact") {
                                          try({
                                                path <- propose_lna(lna_times         = lna_times,
                                                                    lna_pars          = lna_pars,
                                                                    init_start        = stem_object$dynamics$lna_initdist_inds[1],
                                                                    param_update_inds = param_update_inds,
                                                                    stoich_matrix     = stem_object$dynamics$stoich_matrix_lna,
                                                                    forcing_inds      = forcing_inds,
                                                                    forcing_matrix    = forcing_matrix,
                                                                    step_size         = stem_object$dynamics$dynamics_args$step_size,
                                                                    max_attempts      = max_attempts,
                                                                    reject_negatives  = TRUE,
                                                                    lna_pointer       = stem_object$dynamics$lna_pointers$lna_ptr,
                                                                    set_pars_pointer  = stem_object$dynamics$lna_pointers$set_lna_params_ptr)
                                                
                                                
                                          }, silent = TRUE)
                                          
                                          attempt           <- attempt + 1
                                          
                                          if(!is.null(path)) {
                                                census_paths[[k]] <- census_incidence(path$lna_path, census_times, census_interval_inds)
                                                lna_paths[[k]]    <- path$prev_path[prev_inds,]
                                                lna_draws[[k]]    <- path$draws
                                                
                                                # lna_moments, used in testing
                                                # drift_vecs[[k]]   <- path$drift_vecs
                                                # diff_mats[[k]]    <- path$diff_mats
                                          }
                                          
                                    } else if(lna_method == "approx") {
                                          try({
                                                path <- propose_lna_approx(lna_times         = lna_times,
                                                                           lna_pars          = lna_pars,
                                                                           init_start        = stem_object$dynamics$lna_initdist_inds[1],
                                                                           param_update_inds = param_update_inds,
                                                                           stoich_matrix     = stem_object$dynamics$stoich_matrix_lna,
                                                                           forcing_inds      = forcing_inds,
                                                                           forcing_matrix    = forcing_matrix,
                                                                           step_size         = stem_object$dynamics$dynamics_args$step_size,
                                                                           max_attempts      = max_attempts,
                                                                           nsim              = 1,
                                                                           ess_updates       = 1, 
                                                                           ess_warmup        = ess_warmup,
                                                                           lna_pointer       = stem_object$dynamics$lna_pointers$lna_ptr,
                                                                           set_pars_pointer  = stem_object$dynamics$lna_pointers$set_lna_params_ptr)
                                          }, silent = TRUE)
                                          
                                          attempt           <- attempt + 1
                                          
                                          if(!is.null(path)) {
                                                census_paths[[k]] <- census_incidence(path$incid_paths[,,1], census_times, census_interval_inds)
                                                lna_paths[[k]]    <- path$prev_paths[prev_inds,,1]
                                                lna_draws[[k]]    <- path$lna_draws[,,1]
                                                
                                                # lna_moments, used in testing
                                                # drift_vecs[[k]]   <- path$drift_vecs
                                                # diff_mats[[k]]    <- path$diff_mats
                                          }
                                    }
                              }
                              
                              if(is.null(census_paths[[k]])) {
                                    warning("Simulation failed. Increase max attempts per simulation or try different parameter values.")
                              } else {
                                    colnames(census_paths[[k]]) <- c("time", rownames(stem_object$dynamics$flow_matrix_lna))
                                    colnames(lna_paths[[k]]) <- c("time",colnames(stem_object$dynamics$flow_matrix_lna))
                              }
                        }
                        
                        failed_runs  <- which(sapply(lna_paths, is.null))
                        if(length(failed_runs) != 0) {
                              census_paths <- census_paths[-failed_runs]
                              lna_paths    <- lna_paths[-failed_runs]
                              lna_draws    <- lna_draws[-failed_runs]
                        }   

                } else if(method == "ode") {

                        # set the vectors of times when the ODE is evaluated and censused
                        ode_times <- sort(unique(c(t0, census_times, stem_object$dynamics$tcovar[,1], tmax)))

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
                                tcovar_rowinds <- match(round(ode_times, digits = 8),
                                                        round(stem_object$dynamics$tcovar[,1], digits = 8))
                                ode_pars[tcovar_rowinds, tcovar_inds] <- stem_object$dynamics$tcovar[tcovar_rowinds,-1]

                                # zero out forcings if necessary
                                if(!is.null(stem_object$dynamics$dynamics_args$forcings)) {

                                        # get the forcing indices (supplied in the original tcovar matrix)
                                        forcing_inds <- as.logical(match(round(ode_times, digits = 8),
                                                                         round(stem_object$dynamics$dynamics_args$tcovar[,1], digits = 8),
                                                                         nomatch = FALSE))
                                        zero_inds    <- !forcing_inds

                                        # zero out the tcovar elements corresponding to times with no forcings
                                        for(l in seq_along(stem_object$dynamics$dynamics_args$forcings)) {
                                                ode_pars[zero_inds, stem_object$dynamics$dynamics_args$forcings[[l]]$tcovar_name] = 0
                                        }
                                }
                        }
                        
                        # tparam indices and initial values
                        if(!is.null(stem_object$dynamics$tparam)) {
                              
                              # generate indices for time-varying parameters
                              for(s in seq_along(stem_object$dynamics$tparam)) {
                                    stem_object$dynamics$tparam[[s]]$col_ind  <- 
                                          stem_object$dynamics$ode_rates$ode_param_codes[stem_object$dynamics$tparam[[s]]$tparam_name]
                                    stem_object$dynamics$tparam[[s]]$tpar_inds <- 
                                          findInterval(ode_times, tparam[[s]]$times, left.open = F) - 1
                              }
                              
                              # list for saving the time varying parameters for reuse in simulating a dataset in necessary
                              tparam_draws <- vector("list", length = nsim)
                              tparam_times <- sort(unique(unlist(lapply(tparam, function(x) x$times))))
                              tparam_times <- tparam_times[tparam_times >= t0 & tparam_times <= tmax]
                        } else {
                              tparam_draws <- NULL
                              tparam_times <- NULL
                        }

                        # generate some auxilliary objects
                        param_update_inds <- ode_times %in% unique(c(t0, tmax, stem_object$dynamics$tcovar[,1]))
                        census_interval_inds <- findInterval(ode_times, census_times, left.open = T)

                        # objects for storing the ODE paths
                        fixed_parameters <- stem_object$dynamics$fixed_inits & is.null(simulation_parameters)

                        if(!fixed_parameters) {
                                census_paths <- vector(mode = "list", length = nsim)
                                ode_paths    <- vector(mode = "list", length = nsim)
                        } else {
                                census_paths <- vector(mode = "list", length = 1)
                                ode_paths    <- vector(mode = "list", length = 1)
                        }

                        # retrieve the initial state
                        init_state <- ode_pars[1, stem_object$dynamics$ode_initdist_inds + 1]

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

                                                if(!stem_object$dynamics$state_initializer[[s]]$fixed) {

                                                        init_states[,stem_object$dynamics$state_initializer[[s]]$codes] <-
                                                                stem_object$dynamics$strata_sizes[s] * extraDistr::rdirichlet(nsim, stem_object$dynamics$state_initializer[[s]]$prior)

                                                } else {
                                                        init_states[,stem_object$dynamics$state_initializer[[s]]$codes] <-
                                                                matrix(stem_object$dynamics$state_initializer[[s]]$init_states,
                                                                       nrow = nsim,
                                                                       ncol = length(stem_object$dynamics$state_initializer[[s]]$init_states),
                                                                       byrow = TRUE)
                                                }
                                        }
                                }
                        } else {
                                if(stem_object$dynamics$n_strata == 1) {
                                        init_states <- matrix(rep(stem_object$dynamics$state_initializer$init_states, nsim), nrow = nsim, byrow = T)

                                } else {
                                        init_states <- matrix(0.0, nrow = nsim, ncol = length(stem_object$dynamics$comp_codes))
                                        colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                                        for(s in seq_len(stem_object$dynamics$n_strata)) {
                                                init_states[,stem_object$dynamics$state_initializer[[s]]$codes] <-
                                                        matrix(stem_object$dynamics$state_initializer[[s]]$init_states,
                                                               nrow = nsim,
                                                               ncol = length(stem_object$dynamics$state_initializer[[s]]$init_states),
                                                               byrow = TRUE)
                                        }
                                }
                        }

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

                        sim_pars  <- stem_object$dynamics$parameters
                        init_vols <- init_states[1,]

                        for(k in seq_along(census_paths)) {

                                if(!is.null(simulation_parameters)) {
                                        sim_pars <- as.numeric(simulation_parameters[[k]])
                                }

                                if(!stem_object$dynamics$fixed_inits) {
                                        init_vols <- init_states[k,]
                                }

                                # set the parameters and initial volumes
                                pars2lnapars(ode_pars, c(sim_pars, init_vols))
                                
                                # draw values for the time-varying parameters
                                if(!is.null(stem_object$dynamics$tparam)) {
                                      for(s in seq_along(stem_object$dynamics$tparam)) {
                                            
                                            #draw new values
                                            draw_normals(stem_object$dynamics$tparam[[s]]$values)
                                            
                                            # insert the new values into the tcovar matrix
                                            insert_tparam(tcovar    = ode_pars, 
                                                          values    = stem_object$dynamics$tparam[[s]]$draws2par(stem_object$dynamics$parameters,
                                                                                                                 stem_object$dynamics$tparam[[s]]$values),
                                                          col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                                                          tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds)
                                      }
                                      
                                      tparam_draws[[k]] <- lapply(stem_object$dynamics$tparam, function(x) x$values)
                                }

                                path <- NULL

                                try({
                                        path <- integrate_odes(ode_times         = ode_times,
                                                               ode_pars          = ode_pars,
                                                               init_start        = stem_object$dynamics$ode_initdist_inds[1],
                                                               param_update_inds = param_update_inds,
                                                               stoich_matrix     = stem_object$dynamics$stoich_matrix_ode,
                                                               forcing_inds      = forcing_inds,
                                                               forcing_matrix    = forcing_matrix,
                                                               step_size         = stem_object$dynamics$dynamics_args$step_size,
                                                               ode_pointer       = stem_object$dynamics$ode_pointers$ode_ptr,
                                                               set_pars_pointer  = stem_object$dynamics$ode_pointers$set_ode_params_ptr)
                                }, silent = TRUE)

                                if(!is.null(path)) {
                                        census_paths[[k]] <- census_incidence(path$incid_path, census_times, census_interval_inds)
                                        ode_paths[[k]]    <- path$prev_path[match(round(census_times, digits = 8),
                                                                                  round(path$prev_path[,1], digits = 8)),]

                                        colnames(census_paths[[k]]) <- c("time", rownames(stem_object$dynamics$flow_matrix_ode))
                                        colnames(ode_paths[[k]]) <- c("time",colnames(stem_object$dynamics$flow_matrix_ode))
                                }

                                if(is.null(path) & messages) {
                                        warning("Simulation failed. Try different parameter values.")
                                }
                        }

                        failed_runs  <- which(sapply(ode_paths, is.null))
                        if(length(failed_runs) != 0) {
                                census_paths <- census_paths[-failed_runs]
                                ode_paths    <- ode_paths[-failed_runs]
                        }
                }

                if(observations && length(failed_runs)!= nsim) {

                        datasets         <- vector(mode = "list", length = length(census_paths)) # list for storing the datasets
                        measvar_names    <- colnames(stem_object$measurement_process$obsmat)

                        # grab the time-varying covariate values at observation times
                        tcovar_obstimes <- build_census_path(path           = stem_object$dynamics$tcovar,
                                                             census_times   = stem_object$measurement_process$obstimes,
                                                             census_columns = 1:(ncol(stem_object$dynamics$tcovar)-1))
                        colnames(tcovar_obstimes) <- colnames(stem_object$dynamics$tcovar)

                        # if incidence, the incidence codes are not null
                        do_incidence <- !is.null(stem_object$dynamics$incidence_codes) & !(method %in% c("lna", "ode"))

                        # census if computing incidence or if computing prevalence and the obstimes != census_times
                        do_census <- !is.null(census_times) &&
                                !identical(as.numeric(stem_object$measurement_process$obstimes), as.numeric(census_times))

                        # get the indices in the censused matrices for the observation times
                        if(do_census) cens_inds <- findInterval(stem_object$measurement_process$obstimes, census_times)

                        if(method == "gillespie") {

                                # should incidence be computed?
                                if(do_incidence) incidence_codes <- stem_object$dynamics$incidence_codes + 1

                                # column codes in the path matrix for compartments to be censused
                                census_codes    <- c(stem_object$dynamics$comp_codes,
                                                     stem_object$dynamics$incidence_codes) + 2
                                censusmat       <- stem_object$measurement_process$censusmat

                                # initialize simulation parameters
                                sim_pars <- stem_object$dynamics$parameters

                                for(k in seq_len(nsim)) {

                                      # get the state at observation times
                                      if(!is.null(census_paths[[k]])) {
                                            
                                            if(do_census) {
                                                  censusmat[,-1] <- census_paths[[k]][cens_inds, -1]
                                            } else {
                                                  censusmat <- census_paths[[k]]
                                            }
                                            
                                            # get the new simulation parameters if a list was supplied
                                            if(!is.null(simulation_parameters)) sim_pars <- as.numeric(simulation_parameters[[k]])
                                            
                                            # insert the time-varying parameters into the tcovar matrix
                                            if(!is.null(tparam_draws)) {
                                                  
                                                  for(s in seq_along(stem_object$dynamics$tparam)) {
                                                        
                                                        # insert the tparam values into the tcovar matrix
                                                        insert_tparam(tcovar = stem_object$dynamics$tcovar, 
                                                                      values = stem_object$dynamics$tparam[[s]]$draws2par(
                                                                            sim_pars, 
                                                                            tparam_draws[[k]][[s]]),
                                                                      col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                                                                      tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds)
                                                        
                                                        # grab the time-varying covariate values at observation times
                                                        tcovar_obstimes <- build_census_path(path           = stem_object$dynamics$tcovar,
                                                                                             census_times   = stem_object$measurement_process$obstimes,
                                                                                             census_columns = 1:(ncol(stem_object$dynamics$tcovar)-1))
                                                        colnames(tcovar_obstimes) <- colnames(stem_object$dynamics$tcovar)
                                                  }
                                            }
                                            
                                            # simulate the data
                                            datasets[[k]] <- simulate_r_measure(censusmat = censusmat,
                                                                                measproc_indmat = stem_object$measurement_process$measproc_indmat,
                                                                                parameters = sim_pars,
                                                                                constants = stem_object$dynamics$constants,
                                                                                tcovar = tcovar_obstimes,
                                                                                r_measure_ptr = stem_object$measurement_process$meas_pointers$r_measure_ptr)
                                            colnames(datasets[[k]]) <- measvar_names
                                      }
                                }

                        } else if(method == "lna") {

                                # get the objects for simulating from the measurement process
                                measproc_indmat  <- stem_object$measurement_process$measproc_indmat
                                sim_pars         <- stem_object$dynamics$parameters
                                constants        <- stem_object$dynamics$constants
                                tcovar           <- tcovar_obstimes
                                r_measure_ptr    <- stem_object$measurement_process$meas_pointers$r_measure_ptr
                                cens_inds        <- c(0,match(round(stem_object$measurement_process$obstimes, digits = 8),
                                                              round(census_times, digits = 8)) - 1)
                                do_prevalence    <- stem_object$measurement_process$lna_prevalence
                                do_incidence     <- stem_object$measurement_process$lna_incidence
                                obstime_inds     <- stem_object$measurement_process$obstime_inds
                                pathmat          <- stem_object$measurement_process$censusmat
                                flow_matrix_lna  <- stem_object$dynamics$flow_matrix_lna
                                lna_event_inds   <- stem_object$measurement_process$incidence_codes_lna
                                
                                # reinitialize the tparam indices if necessary
                                if(!is.null(stem_object$dynamics$tparam)) {
                                      # reinitialize the tparam indices if necessary
                                      for(s in seq_along(stem_object$dynamics$tparam)) {
                                            stem_object$dynamics$tparam[[s]]$col_ind  <- 
                                                  stem_object$dynamics$tcovar_codes[stem_object$dynamics$tparam[[s]]$tparam_name]
                                            stem_object$dynamics$tparam[[s]]$tpar_inds <- 
                                                  findInterval(stem_object$dynamics$tcovar[,1], tparam[[s]]$times, left.open = F) - 1
                                      }
                                }

                                for(k in seq_along(census_paths)) {

                                        # fill out the census matrix
                                        census_lna(path                = census_paths[[k]],
                                                   census_path         = pathmat,
                                                   census_inds         = cens_inds,
                                                   lna_event_inds      = lna_event_inds,
                                                   flow_matrix_lna     = flow_matrix_lna,
                                                   do_prevalence       = do_prevalence,
                                                   init_state          = init_states[k,],
                                                   forcing_matrix      = forcing_matrix)

                                        if(!is.null(simulation_parameters)) sim_pars <- as.numeric(simulation_parameters[[k]])
                                        
                                        # insert the time-varying parameters into the tcovar matrix
                                        if(!is.null(tparam_draws)) {
                                              
                                              for(s in seq_along(stem_object$dynamics$tparam)) {
                                                    
                                                    # insert the tparam values into the tcovar matrix
                                                    insert_tparam(tcovar = stem_object$dynamics$tcovar, 
                                                                  values = stem_object$dynamics$tparam[[s]]$draws2par(
                                                                        sim_pars, 
                                                                        tparam_draws[[k]][[s]]),
                                                                  col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                                                                  tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds)
                                                    
                                                    # grab the time-varying covariate values at observation times
                                                    tcovar_obstimes <- build_census_path(path           = stem_object$dynamics$tcovar,
                                                                                         census_times   = stem_object$measurement_process$obstimes,
                                                                                         census_columns = 1:(ncol(stem_object$dynamics$tcovar)-1))
                                                    colnames(tcovar_obstimes) <- colnames(stem_object$dynamics$tcovar)
                                              }
                                        }

                                        # simulate the dataset
                                        datasets[[k]] <- simulate_r_measure(pathmat,
                                                                            measproc_indmat,
                                                                            sim_pars,
                                                                            constants,
                                                                            tcovar_obstimes,
                                                                            r_measure_ptr)
                                        colnames(datasets[[k]]) <- measvar_names
                                }

                        } else if(method == "ode") {

                                # get the objects for simulating from the measurement process
                                measproc_indmat  <- stem_object$measurement_process$measproc_indmat
                                sim_pars         <- stem_object$dynamics$parameters
                                constants        <- stem_object$dynamics$constants
                                tcovar           <- tcovar_obstimes
                                r_measure_ptr    <- stem_object$measurement_process$meas_pointers$r_measure_ptr
                                cens_inds        <- c(0,match(round(stem_object$measurement_process$obstimes, digits = 8),
                                                              round(census_times, digits = 8)) - 1)
                                do_prevalence    <- stem_object$measurement_process$ode_prevalence
                                do_incidence     <- stem_object$measurement_process$ode_incidence
                                obstime_inds     <- stem_object$measurement_process$obstime_inds
                                pathmat          <- stem_object$measurement_process$censusmat
                                flow_matrix_ode  <- stem_object$dynamics$flow_matrix_ode
                                ode_event_inds   <- stem_object$measurement_process$incidence_codes_ode

                                if(!is.null(stem_object$dynamics$tparam)) {
                                      # reinitialize the tparam indices if necessary
                                      for(s in seq_along(stem_object$dynamics$tparam)) {
                                            stem_object$dynamics$tparam[[s]]$col_ind  <- 
                                                  stem_object$dynamics$tcovar_codes[stem_object$dynamics$tparam[[s]]$tparam_name]
                                            stem_object$dynamics$tparam[[s]]$tpar_inds <- 
                                                  findInterval(stem_object$dynamics$tcovar[,1], tparam[[s]]$times, left.open = F) - 1
                                      }
                                }
                                
                                if(!fixed_parameters) {

                                        for(k in seq_along(census_paths)) {

                                                # fill out the census matrix
                                                census_lna(path                = census_paths[[k]],
                                                           census_path         = pathmat,
                                                           census_inds         = cens_inds,
                                                           lna_event_inds      = ode_event_inds,
                                                           flow_matrix_lna     = flow_matrix_ode,
                                                           do_prevalence       = do_prevalence,
                                                           init_state          = init_states[k,],
                                                           forcing_matrix      = forcing_matrix)

                                                if(!is.null(simulation_parameters)) sim_pars <- as.numeric(simulation_parameters[[k]])

                                                # insert the time-varying parameters into the tcovar matrix
                                                if(!is.null(tparam_draws)) {
                                                      
                                                      for(s in seq_along(stem_object$dynamics$tparam)) {
                                                            
                                                            # insert the tparam values into the tcovar matrix
                                                            insert_tparam(tcovar = stem_object$dynamics$tcovar, 
                                                                          values = stem_object$dynamics$tparam[[s]]$draws2par(
                                                                                sim_pars, 
                                                                                tparam_draws[[k]][[s]]),
                                                                          col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                                                                          tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds)
                                                            
                                                            # grab the time-varying covariate values at observation times
                                                            tcovar_obstimes <- build_census_path(path           = stem_object$dynamics$tcovar,
                                                                                                 census_times   = stem_object$measurement_process$obstimes,
                                                                                                 census_columns = 1:(ncol(stem_object$dynamics$tcovar)-1))
                                                            colnames(tcovar_obstimes) <- colnames(stem_object$dynamics$tcovar)
                                                      }
                                                }
                                                
                                                # simulate the dataset
                                                datasets[[k]] <- simulate_r_measure(pathmat,
                                                                                    measproc_indmat,
                                                                                    sim_pars,
                                                                                    constants,
                                                                                    tcovar_obstimes,
                                                                                    r_measure_ptr)
                                                colnames(datasets[[k]]) <- measvar_names
                                        }
                                } else {
                                        # fill out the census matrix
                                        census_lna(path                = census_paths[[1]],
                                                   census_path         = pathmat,
                                                   census_inds         = cens_inds,
                                                   lna_event_inds      = ode_event_inds,
                                                   flow_matrix_lna     = flow_matrix_ode,
                                                   do_prevalence       = do_prevalence,
                                                   init_state          = init_states[1,],
                                                   forcing_matrix      = forcing_matrix)

                                        for(k in seq_len(nsim)) {

                                              # insert the time-varying parameters into the tcovar matrix
                                              if(!is.null(tparam_draws)) {
                                                    
                                                    for(s in seq_along(stem_object$dynamics$tparam)) {
                                                          
                                                          # insert the tparam values into the tcovar matrix
                                                          insert_tparam(tcovar = stem_object$dynamics$tcovar, 
                                                                        values = stem_object$dynamics$tparam[[s]]$draws2par(
                                                                              sim_pars, 
                                                                              tparam_draws[[k]][[s]]),
                                                                        col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                                                                        tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds)
                                                          
                                                          # grab the time-varying covariate values at observation times
                                                          tcovar_obstimes <- build_census_path(path           = stem_object$dynamics$tcovar,
                                                                                               census_times   = stem_object$measurement_process$obstimes,
                                                                                               census_columns = 1:(ncol(stem_object$dynamics$tcovar)-1))
                                                          colnames(tcovar_obstimes) <- colnames(stem_object$dynamics$tcovar)
                                                    }
                                              }
                                              
                                                # simulate the dataset
                                                datasets[[k]] <- simulate_r_measure(pathmat,
                                                                                    measproc_indmat,
                                                                                    sim_pars,
                                                                                    constants,
                                                                                    tcovar_obstimes,
                                                                                    r_measure_ptr)
                                                colnames(datasets[[k]]) <- measvar_names
                                        }
                                }
                        }
                } else if(observations && length(failed_runs == nsim)) {
                      datasets = list(NULL)
                }

                stem_simulations <- list(paths = NULL, datasets = NULL, subject_paths = NULL)

                if(paths & is.null(census_times)) {
                        stem_simulations$paths <- paths_full
                } else if(paths & !is.null(census_times)) {
                        stem_simulations$paths <- census_paths
                        if(method == "lna") stem_simulations$natural_paths <- lna_paths
                        if(method == "ode") stem_simulations$natural_paths <- ode_paths
                } else {
                        stem_simulations$paths <- NULL
                }
                if(observations)  stem_simulations$datasets      <- datasets
                if(subject_paths) stem_simulations$subject_paths <- subject_paths
                if(method == "lna") stem_simulations$lna_draws   <- lna_draws
                stem_simulations$failed_runs <- failed_runs

                # lna_moments used in testing
                #                 if(method == "lna") stem_simulations$lna_moments <- list(drift_vecs = drift_vecs,
                #                                                                          diff_mats = diff_mats)

                return(stem_simulations)
        }