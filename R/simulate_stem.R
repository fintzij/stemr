#' Simulations from a stochastic epidemic model.
#'
#' @param nsim number of realizations to simulate
#' @param simulation_parameters optional list of vectors of simulation
#'   parameters. If NULL, paths are simulated using the parameters specified in
#'   the stem object.
#' @param lna_draws optional list of matrices of lna draws
#' @param tparam_draws optional list of lists of time-varying parameter draws
#' @param paths Should population-level paths at census times be returned?
#' @param full_paths Should complete population level paths be returned (only
#'   for Gillespie), defaults to FALSE
#' @param observations Should simulated observations be returned? Requires that
#'   a measurement process be defined in the stem object.
#' @param method either "gillespie" if simulating via Gillespie's direct method,
#'   "lna" if simulating paths via the linear noise approximation, or "ode" if
#'   simulating paths of the deterministic limit of the underlying Markov jump
#'   process.
#' @param tmax the time at which simulation of the system is terminated.
#'   Defaults to the last observation time if not supplied.
#' @param census_times vector of times at which compartment counts should be
#'   recorded.
#' @param max_attempts maximum number of times to attempt simulating an LNA path
#'   before aborting due to the moments of the path being degenerate at some
#'   point
#' @param lna_method defaults to "exact". If "approx", an initial path is
#'   proposed by resampling perturbations that lead to negative increments or
#'   volumes. The initial path is then updated via elliptical slice sampling
#' @param ess_warmup number of elliptical slice sampling updates before the lna
#'   sample is saved
#' @param messages should a message be printed when parsing the rates?
#' @param stem_object stem object list
#' @param lna_bracket_width initial elliptical slice sampling bracket width to
#'   be used if lna_method == "approx"
#'
#' @return Returns a list with the simulated paths, subject-level paths, and/or
#'   datasets. If \code{paths = FALSE} and \code{observations = FALSE}, or if
#'   \code{paths = TRUE} and \code{observations = TRUE}, a list \code{nsim stem}
#'   paths and datasets, each returned either as a list or array, is returned.
#'
#'   If \code{paths = TRUE} and \code{observations = FALSE}, a list or array of
#'   simulated population-level paths is returned.
#'
#'   If \code{paths = FALSE} and \code{observations = TRUE}, a list or array of
#'   simulated datasets is returned.
#' @export
simulate_stem <-
    function(stem_object,
             nsim = 1,
             simulation_parameters = NULL,
             lna_draws = NULL,
             tparam_draws = NULL,
             tparam_values = NULL,
             paths = FALSE,
             full_paths = FALSE,
             observations = FALSE,
             method = "gillespie",
             tmax = NULL,
             census_times = NULL,
             max_attempts = 500,
             lna_method = "exact",
             lna_bracket_width = 2 * pi,
             ess_warmup = 100,
             messages = TRUE) {

        # ensure that the method is correctly specified
        if (!method %in% c("gillespie", "lna", "ode")) {
            stop("The simulation method must either be 'gillespie', 'lna', or 'ode'.")
        }

        # make sure the object was appropriately compiled
        if (method == "gillespie" &
            is.null(stem_object$dynamics$rate_ptrs)) {
            stop("Exact rates not compiled.")
        } else if (method == "lna" &
                   is.null(stem_object$dynamics$lna_pointers)) {
            stop("LNA not compiled.")
        } else if (method == "ode" &
                   is.null(stem_object$dynamics$ode_pointers)) {
            stop("ODE not compiled.")
        }

        # check that the stem_dynamics are supplied
        if (is.null(stem_object$dynamics)) {
            stop("The stochastic epidemic model dynamics must be specified.")
        }

        # check that a measurement process is supplied if one is required
        if (is.null(stem_object$measurement_process) &
            observations) {
            stop(
                "In order to simulate a dataset, a measurement process must be supplied in the stem_object."
            )
        }

        # check that if simulation parameters are supplied, they are supplied in a list of named vectors of
        # length equal to the number of simulations. Also check that the parameters are supplied in the same
        # order as the parameters in the stem object.
        if (!is.null(simulation_parameters)) {
            if (!is.list(simulation_parameters)) {
                stop("Simulation parameters must be supplied as a list of vectors.")
            }

            if (length(simulation_parameters) != nsim) {
                stop(
                    "If simulation parameters are supplied, there must be an equal number of vectors in the list as requested paths to be simulated."
                )
            }

            if (!all(sapply(simulation_parameters, function(x)
                all(names(x) == 
                    names(stem_object$dynamics$parameters)[!grepl("_0", names(stem_object$dynamics$parameters))]
                )))) {
                stop(
                    "The simulation parameters list must consist of named vectors with elements given in the same order as the parameters in the stem object."
                )
            }
        }

        if (!is.null(tparam_draws) & !is.null(tparam_values)) {
            warning(
                "Draws and values for time varying parameters were both specified. Values will take precedence."
            )
        }

        if (!is.null(tparam_draws) & length(tparam_draws) != nsim) {
            stop(
                "The number of time varying parameter draw lists must be equal to the number of simulations."
            )
        }

        if (!is.null(tparam_values) &
            length(tparam_values) != nsim) {
            stop(
                "The number of time varying parameter value lists must be equal to the number of simulations."
            )
        }

        if (paths == observations) {
            paths <- TRUE
            observations <-
                ifelse(is.null(stem_object$measurement_process),
                       FALSE,
                       TRUE)
        }

        # t0 and tmax if not supplied
        t0                     <- stem_object$dynamics$t0
        if (is.null(tmax))
            tmax <- stem_object$dynamics$tmax

        if (is.null(stem_object$dynamics$timestep)) {
            timestep <- 1
        } else {
            timestep <- stem_object$dynamics$timestep
        }

        # make sure that there exists a vector of census times
        if (is.null(census_times)) {
            if (is.null(stem_object$measurement_process$obstimes)) {
                census_times <-
                    as.numeric(unique(c(
                        t0, seq(t0, tmax, by = timestep), tmax
                    )))

            } else {
                census_times <-
                    as.numeric(unique(
                        c(
                            t0,
                            stem_object$measurement_process$obstimes,
                            tmax
                        )
                    ))
            }

        } else {
            census_times <- as.numeric(sort(unique(c(
                t0, census_times, tmax
            ))))
        }

        if (!is.numeric(stem_object$dynamics$parameters)) {
            stop("The model parameters must be a named numeric vector.")
        }

        if (!is.numeric(stem_object$dynamics$initdist_params)) {
            stop("The initial distribution parameters must be a named numeric vector.")
        }

        if (!is.numeric(stem_object$dynamics$constants)) {
            stop("The model constants must be a named numeric vector.")
        }

        # build the time varying covariate matrix (includes, at a minimum, the endpoints of the simulation interval)
        # if timestep is null, there are no time-varying covariates
        if (method == "gillespie") {

            # if any of t0, tmax, or a timestep was supplied,
            # check if they differ from the parameters supplied in the stem_object$dynamics.
            # if they differ, reconstruct the tcovar matrix and associated objects
            rebuild_tcovar <- (t0 != stem_object$dynamics$t0) ||
                (tmax != stem_object$dynamics$tcovar[nrow(stem_object$dynamics$tcovar), 1]) ||
                timestep != stem_object$dynamics$timestep ||
                !identical(stem_object$dynamics$tcovar[, 1], census_times)

            if (rebuild_tcovar) {
                # check for explicit time reference
                no_TIME_ref <-
                    !any(sapply(stem_object$dynamics$rates, function(x)
                        grepl("TIME", x[["unparsed"]])))

                # rebuild the time-varying covariate matrix so that it contains the census intervals
                stem_object$dynamics$tcovar <-
                    build_tcovar_matrix(
                        tcovar       = stem_object$dynamics$dynamics_args$tcovar,
                        tparam       = stem_object$dynamics$tparam,
                        forcings     = stem_object$dynamics$forcings,
                        parameters   = stem_object$dynamics$parameters,
                        timestep     = timestep,
                        census_times = census_times,
                        t0           = t0,
                        tmax         = tmax,
                        messages     = messages
                    )

                stem_object$dynamics$tcovar_codes        <-
                    seq_len(ncol(stem_object$dynamics$tcovar) - 1)
                names(stem_object$dynamics$tcovar_codes) <-
                    colnames(stem_object$dynamics$tcovar)[2:ncol(stem_object$dynamics$tcovar)]
                stem_object$dynamics$n_tcovar            <-
                    ncol(stem_object$dynamics$tcovar) - 1
                stem_object$dynamics$tcovar_changemat    <-
                    build_tcovar_changemat(
                        tcovar = stem_object$dynamics$tcovar,
                        tparam = stem_object$dynamics$tparam,
                        forcings = stem_object$dynamics$forcings,
                        no_TIME_ref = no_TIME_ref
                    )
                stem_object$dynamics$tcovar_adjmat       <-
                    build_tcovar_adjmat(
                        rates        = stem_object$dynamics$rates,
                        tcovar_codes = stem_object$dynamics$tcovar_codes,
                        forcings     = stem_object$dynamics$forcings
                    )

                # zero out forcings if necessary
                if (!is.null(stem_object$dynamics$dynamics_args$forcings)) {
                    # get the forcing indices (supplied in the original tcovar matrix)
                    forcing_inds <-
                        vector("logical",
                               nrow(stem_object$dynamics$tcovar))
                    for (f in seq_along(stem_object$dynamics$forcings)) {
                        forcing_inds <-
                            forcing_inds |
                            stem_object$dynamics$tcovar[, stem_object$dynamics$forcings[[f]]$tcovar_name] != 0
                    }

                } else {
                    forcing_inds   <- rep(FALSE, nrow(stem_object$dynamics$tcovar))
                }

            } else {
                # Get the forcing indices if there are forcings
                if (!is.null(stem_object$dynamics$dynamics_args$forcings)) {
                    # get the forcing indices (supplied in the original tcovar matrix)
                    forcing_inds <-
                        vector("logical",
                               nrow(stem_object$dynamics$tcovar))
                    for (f in seq_along(stem_object$dynamics$forcings)) {
                        forcing_inds <-
                            forcing_inds |
                            stem_object$dynamics$tcovar[, stem_object$dynamics$forcings[[f]]$tcovar_name] != 0
                    }

                } else {
                    forcing_inds   <- rep(FALSE, nrow(stem_object$dynamics$tcovar))
                }
            }

            # generate or copy the initial states
            if (stem_object$dynamics$fixed_inits) {
                # if all initial states are fixed, just copy the initial compartment counts
                init_states <- matrix(rep(
                    as.numeric(stem_object$dynamics$initdist_params),
                    nsim
                ),
                nrow = nsim,
                byrow = TRUE)
                colnames(init_states) <-
                    names(stem_object$dynamics$comp_codes)

            } else {
                if (stem_object$dynamics$n_strata == 1) {
                    # simulate the initial compartment counts
                    init_states <-
                        t(as.matrix(
                            rmultinom(
                                nsim,
                                stem_object$dynamics$popsize,
                                stem_object$dynamics$initdist_priors
                            )
                        ))
                    colnames(init_states) <-
                        names(stem_object$dynamics$comp_codes)

                } else if (stem_object$dynamics$n_strata > 1) {
                    # generate the matrix of initial compartment counts
                    init_states <-
                        matrix(
                            0,
                            nrow = nsim,
                            ncol = length(stem_object$dynamics$comp_codes)
                        )
                    colnames(init_states) <-
                        names(stem_object$dynamics$comp_codes)

                    for (s in seq_len(stem_object$dynamics$n_strata)) {
                        if (!stem_object$dynamics$initializer[[s]]$fixed) {
                            if (stem_object$dynamics$iniitializer[[s]]$dist == "multinom") {
                                init_states[, stem_object$dynamics$initializer[[s]]$codes] <-
                                    as.matrix(t(
                                        rmultinom(
                                            nsim,
                                            stem_object$dynamics$strata_sizes[s],
                                            stem_object$dynamics$initializer[[s]]$prior
                                        )
                                    ))
                            } else {
                                init_states[, stem_object$dynamics$initializer[[s]]$codes] <-
                                    extraDistr::rdirmnom(
                                        nsim,
                                        stem_object$dynamics$strata_sizes[s],
                                        ifelse(
                                            stem_object$dynamics$initializer[[s]]$prior != 0,
                                            stem_object$dynamics$initializer[[s]]$prior,
                                            .Machine$double.eps
                                        )
                                    )
                            }

                        } else {
                            init_states[, stem_object$dynamics$initializer[[s]]$codes] <-
                                matrix(
                                    as.numeric(
                                        stem_object$dynamics$initializer[[s]]$init_states
                                    ),
                                    nrow = nsim,
                                    ncol = length(
                                        stem_object$dynamics$initializer[[s]]$init_states
                                    ),
                                    byrow = T
                                )
                        }
                    }
                }
            }

            # if there are artificial incidence compartments, copy the
            # incidence counts and add them to the initial state matrix
            if (!is.null(stem_object$dynamics$incidence_codes)) {
                # init_incid <- init_states[, stem_object$dynamics$incidence_sources + 1, drop = FALSE]
                init_incid <-
                    matrix(
                        0,
                        nrow = nrow(init_states),
                        ncol = length(stem_object$dynamics$incidence_codes)
                    )
                colnames(init_incid) <-
                    names(stem_object$dynamics$incidence_codes)
                init_states <- cbind(init_states, init_incid)
            }

            # guess the initial dimensions. need an extra column for event times and another for event IDs.
            if (stem_object$dynamics$progressive & any(stem_object$dynamics$absorbing_states)) {
                if (stem_object$dynamics$n_strata == 1) {
                    init_dims <-
                        c(
                            n_rows = sum(
                                stem_object$dynamics$popsize * stem_object$dynamics$n_compartments
                            ),
                            n_cols = ncol(stem_object$dynamics$flow_matrix) + 2
                        )
                } else {
                    init_dims <-
                        c(
                            n_rows = sum(
                                stem_object$dynamics$strata_sizes * sapply(
                                    sapply(stem_object$dynamics$initializer, "[[", 4),
                                    length
                                )
                            ),
                            n_cols = ncol(stem_object$dynamics$flow_matrix) + 2
                        )

                }
            } else {
                if (stem_object$dynamics$n_strata == 1) {
                    init_dims <-
                        c(
                            n_rows = sum(
                                stem_object$dynamics$popsize * stem_object$dynamics$n_compartments
                            ) * 3,
                            n_cols = ncol(stem_object$dynamics$flow_matrix) + 2
                        )
                } else if (stem_object$dynamics$n_strata > 1) {
                    init_dims <-
                        c(
                            n_rows = sum(
                                stem_object$dynamics$strata_sizes * sapply(
                                    sapply(stem_object$dynamics$initializer, "[[", 4),
                                    length
                                )
                            ) * 3,
                            n_cols = ncol(stem_object$dynamics$flow_matrix) + 2
                        )
                }
            }

            # not too big
            init_dims[1] <- pmin(init_dims[1], 2 ^ 23)

            # make the initial dimensions a little bigger (round up to nearest power of 2)
            p <- 1
            while (2 ^ p < init_dims[1] & p < 1e6) {
                p <- p + 1
                if (2 ^ p > init_dims[1])
                    init_dims[1] <- 2 ^ p
            }

            # get the compartment names
            path_colnames <-
                c("time", "event", c(
                    names(stem_object$dynamics$comp_codes),
                    names(stem_object$dynamics$incidence_codes)
                ))

            # generate forcing matrix
            if (!is.null(stem_object$dynamics$dynamics_args$forcings)) {
                forcings <- stem_object$dynamic_args$forcings

                # names and indices
                forcing_tcovars   <-
                    sapply(forcings, function(x)
                        x$tcovar_name)
                forcing_tcov_inds <-
                    match(forcing_tcovars,
                          colnames(stem_object$dynamics$tcovar)) - 1
                forcing_events    <-
                    c(sapply(forcings, function(x)
                        paste0(x$from, "2", x$to)))

                # matrix indicating which compartments are involved in which forcings in and out
                forcings_out <- matrix(
                    0.0,
                    nrow = length(path_colnames) - 2,
                    ncol = length(forcings),
                    dimnames = list(path_colnames[-c(1:2)], forcing_tcovars)
                )

                forcing_transfers <- array(
                    0.0,
                    dim = c(
                        length(path_colnames) - 2,
                        length(path_colnames) - 2,
                        length(forcings)
                    ),
                    dimnames = list(path_colnames[-c(1:2)], path_colnames[-c(1:2)], forcing_tcovars)
                )

                for (s in seq_along(forcings)) {
                    forcings_out[forcings[[s]]$from, s] <- 1

                    for (t in seq_along(forcings[[s]]$from)) {
                        forcing_transfers[forcings[[s]]$from[t], forcings[[s]]$from[t], s] <-
                            -1
                        forcing_transfers[forcings[[s]]$to[t], forcings[[s]]$from[t], s]    <-
                            1
                    }
                }

            } else {
                forcing_tcovars   <- character(0L)
                forcing_tcov_inds <- integer(0L)
                forcing_events    <- character(0L)
                forcings_out      <- matrix(0.0, nrow = 0, ncol = 0)
                forcing_transfers <- array(0.0, dim = c(0, 0, 0))
            }

            # vector of model parameters
            sim_pars <- stem_object$dynamics$parameters

            # tparam indices and initial values
            if (!is.null(stem_object$dynamics$tparam)) {
                # generate indices for time-varying parameters
                for (s in seq_along(stem_object$dynamics$tparam)) {
                    stem_object$dynamics$tparam[[s]]$col_ind  <-
                        stem_object$dynamics$tcovar_codes[stem_object$dynamics$tparam[[s]]$tparam_name]
                    stem_object$dynamics$tparam[[s]]$tpar_inds <-
                        findInterval(
                            stem_object$dynamics$tcovar[, 1],
                            stem_object$dynamics$tparam[[s]]$times,
                            left.open = F
                        ) - 1
                    stem_object$dynamics$tparam[[s]]$tpar_inds[stem_object$dynamics$tparam[[s]]$tpar_inds == -1] <-
                        0
                }

                # list for saving the time varying parameters for reuse in simulating a dataset in necessary
                tparam_times <-
                    sort(unique(unlist(
                        lapply(stem_object$dynamics$tparam, function(x)
                            x$times)
                    )))
                tparam_times <-
                    tparam_times[tparam_times >= t0 & tparam_times <= tmax]

                if (is.null(tparam_values) &
                    is.null(tparam_draws)) {
                    tpar_list <-
                        list(lapply(stem_object$dynamics$tparam, function(x)
                            numeric(length(x$values))))

                    tparam_draws  <- rep(tpar_list, nsim)
                    tparam_values <- rep(tpar_list, nsim)

                    # compute the tparam values
                    for (n in seq_len(nsim)) {
                        for (m in seq_along(stem_object$dynamics$tparam)) {
                            # grab parameters
                            if (!is.null(simulation_parameters)) {
                                sim_pars <- simulation_parameters[[n]]
                                class(sim_pars) <- "numeric"
                            }

                            # draw values
                            draw_normals(tparam_draws[[n]][[m]])

                            # compute values
                            tparam_values[[n]][[m]] <-
                                stem_object$dynamics$tparam[[m]]$draws2par(sim_pars,
                                                                           tparam_draws[[n]][[m]])
                        }
                    }

                } else if (is.null(tparam_values) &
                           !is.null(tparam_draws)) {
                    tpar_list <-
                        list(lapply(stem_object$dynamics$tparam, function(x)
                            numeric(length(x$values))))

                    tparam_values <- rep(tpar_list, nsim)

                    # compute the tparam values
                    for (n in seq_len(nsim)) {
                        for (m in seq_along(stem_object$dynamics$tparam)) {
                            # grab parameters
                            if (!is.null(simulation_parameters)) {
                                sim_pars <- simulation_parameters[[n]]
                                class(sim_pars) <- "numeric"
                            }

                            # compute values
                            tparam_values[[n]][[m]] <-
                                stem_object$dynamics$tparam[[m]]$draws2par(sim_pars,
                                                                           tparam_draws[[n]][[m]])
                        }
                    }
                }

            } else {
                tparam_draws <- NULL
                tparam_times <- NULL
            }

            # initialize the list of paths
            if (full_paths) paths_full <- vector(mode = "list", length = nsim)

            census_paths    <- vector(mode = "list", length = nsim)
            census_colnames <-
                c("time", c(
                    names(stem_object$dynamics$comp_codes),
                    names(stem_object$dynamics$incidence_codes)
                ))

            # add 2 to the codes b/c 'time' and 'event' are in the full path
            census_codes  <-
                c(stem_object$dynamics$comp_codes,
                  stem_object$dynamics$incidence_codes) + 2
            get_incidence <-
                !is.null(stem_object$dynamics$incidence_codes)

            if (get_incidence) {
                incidence_codes       <- stem_object$dynamics$incidence_codes + 1
                census_incidence_rows <-
                    rep(list(seq_along(census_times) - 1),
                        length(stem_object$dynamics$incidence_codes))
            }

            for (k in seq_len(nsim)) {

                attempt <- 0
                path_full <- NULL

                if (!is.null(simulation_parameters)) {
                    sim_pars <- simulation_parameters[[k]]
                    class(sim_pars) <- "numeric"
                }

                # draw new time-varying parameters if necessary
                if (!is.null(stem_object$dynamics$tparam)) {
                    for (s in seq_along(stem_object$dynamics$tparam)) {
                        # insert the new values into the tcovar matrix
                        insert_tparam(
                            tcovar = stem_object$dynamics$tcovar,
                            values = tparam_values[[k]][[s]],
                            col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                            tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds
                        )

                    }
                }

                while (is.null(path_full) &
                       attempt < max_attempts) {
                    try({
                        path_full <-
                            simulate_gillespie(
                                flow              = stem_object$dynamics$flow_matrix,
                                parameters        = sim_pars,
                                constants         = stem_object$dynamics$constants,
                                tcovar            = stem_object$dynamics$tcovar,
                                t_max             = max(census_times),
                                init_states       = init_states[k,],
                                rate_adjmat       = stem_object$dynamics$rate_adjmat,
                                tcovar_adjmat     = stem_object$dynamics$tcovar_adjmat,
                                tcovar_changemat  = stem_object$dynamics$tcovar_changemat,
                                init_dims         = init_dims,
                                forcing_inds      = forcing_inds,
                                forcing_tcov_inds = forcing_tcov_inds,
                                forcings_out      = forcings_out,
                                forcing_transfers = forcing_transfers,
                                rate_ptr          = stem_object$dynamics$rate_ptrs[[1]]
                            )
                    }, silent = TRUE)
                    attempt <- attempt + 1
                }

                # get the census path
                if (!is.null(path_full)) {
                    census_paths[[k]] <- build_census_path(
                        path = path_full,
                        census_times = census_times,
                        census_columns = census_codes
                    )

                    # compute incidence if required. n.b. add 1 to the incidence codes b/c 'time' is in the census path
                    if (get_incidence)
                        compute_incidence(
                            censusmat = census_paths[[k]],
                            col_inds  = incidence_codes,
                            row_inds  = census_incidence_rows
                        )

                    # assign column names
                    colnames(census_paths[[k]]) <- census_colnames

                    # save the full path if desired
                    if (full_paths) {
                        paths_full[[k]] <-
                            path_full[, c(1:(
                                2 + length(stem_object$dynamics$comp_codes)
                            ))]
                        colnames(paths_full[[k]]) <-
                            c("time",
                              "event_code",
                              names(stem_object$dynamics$comp_codes))
                    }
                }
            }

            failed_runs  <- which(sapply(census_paths, is.null))
            if (length(failed_runs) != 0) {
                census_paths <- census_paths[-failed_runs]
                if (full_paths)
                    paths_full <- paths_full[-failed_runs]
            }

        } else if (method == "lna") {
            # pull out logical vector for which rates are of order > 1
            higher_order_rates <-
                sapply(stem_object$dynamics$rates, function(x)
                    x$higher_order)

            # set the vectors of times when the LNA is evaluated and censused
            lna_times <- sort(unique(
                c(
                    t0,
                    census_times,
                    seq(t0, tmax, by = stem_object$dynamics$timestep),
                    stem_object$dynamics$tcovar[, 1],
                    tmax
                )
            ))

            lna_census_times <-
                lna_times[lna_times >= t0 & lna_times <= tmax]

            # generate the matrix of parameters, constants, and time-varying covariates
            lna_pars  <- matrix(
                0.0,
                nrow = length(lna_times),
                ncol = length(stem_object$dynamics$lna_rates$lna_param_codes),
                dimnames = list(
                    NULL,
                    names(stem_object$dynamics$lna_rates$lna_param_codes)
                )
            )

            # insert parameters, constants, and time-varying covariates
            parameter_inds <-
                setdiff(
                    stem_object$dynamics$param_codes,
                    stem_object$dynamics$lna_initdist_inds
                )
            constant_inds  <-
                length(stem_object$dynamics$param_codes) + seq_along(stem_object$dynamics$const_codes) - 1
            tcovar_inds    <-
                length(stem_object$dynamics$param_codes) + length(constant_inds) + seq_along(stem_object$dynamics$tcovar_codes) - 1

            lna_pars[, parameter_inds + 1] <- matrix(0.0,
                                                     nrow = nrow(lna_pars),
                                                     ncol = length(parameter_inds))
            lna_pars[, constant_inds + 1]  <- matrix(0.0,
                                                     nrow = nrow(lna_pars),
                                                     ncol = length(constant_inds))

            pars2lnapars2(
                lnapars = lna_pars,
                parameters = c(as.numeric(
                    stem_object$dynamics$parameters[parameter_inds + 1]
                )),
                c_start = 0
            )

            pars2lnapars2(
                lnapars = lna_pars,
                parameters = as.numeric(stem_object$dynamics$constants),
                c_start = constant_inds[1]
            )

            # Generate forcing indices and update the LNA parameter matrix with forcings
            forcing_inds <- rep(FALSE, length(lna_census_times))

            if (!is.null(stem_object$dynamics$dynamics_args$tcovar)) {
                tcovar_rowinds <-
                    findInterval(lna_times,
                                 stem_object$dynamics$tcovar[, 1],
                                 left.open = F)
                lna_pars[tcovar_rowinds, tcovar_inds + 1] <-
                    stem_object$dynamics$tcovar[tcovar_rowinds, -1]

                # zero out forcings if necessary
                if (!is.null(stem_object$dynamics$dynamics_args$forcings)) {
                    # get the forcing indices (supplied in the original tcovar matrix)
                    for (f in seq_along(stem_object$dynamics$forcings)) {
                        forcing_inds <-
                            forcing_inds |
                            stem_object$dynamics$tcovar[stem_object$dynamics$tcovar[, 1] %in% lna_census_times,
                                                        stem_object$dynamics$forcings[[f]]$tcovar_name] != 0
                    }

                    zero_inds    <- !forcing_inds

                    # zero out the tcovar elements corresponding to times with no forcings
                    # zero out the tcovar elements corresponding to times with no forcings
                    for (l in seq_along(stem_object$dynamics$dynamics_args$forcings)) {
                        lna_pars[zero_inds, stem_object$dynamics$dynamics_args$forcings[[l]]$tcovar_name]  = 0
                    }
                }
            }

            # grab parameters
            sim_pars  <-
                stem_object$dynamics$parameters[parameter_inds + 1]

            # simulate the LNA paths
            census_paths <- vector(mode = "list", length = nsim)
            lna_paths    <- vector(mode = "list", length = nsim)
            if (is.null(lna_draws)) {
                lna_draws <-
                    lapply(seq_len(nsim),
                           function(x)
                               matrix(
                                   rnorm(
                                       ncol(stem_object$dynamics$stoich_matrix_lna) * (length(lna_times) - 1)
                                   ),
                                   nrow = ncol(stem_object$dynamics$stoich_matrix_lna)
                               ))
            }

            # retrieve the initial state
            init_state <-
                lna_pars[1, stem_object$dynamics$lna_initdist_inds + 1]

            # object for sampling initial states
            n_strata <- stem_object$dynamics$n_strata
            initializer <- stem_object$dynamics$initializer
            initdist_objects <-
                vector("list", length = stem_object$dynamics$n_strata)

            if (n_strata == 1) {
                comp_size_vec <- stem_object$dynamics$constants["popsize"]
            } else {
                comp_size_vec <-
                    stem_object$dynamics$constants[paste0("popsize_", sapply(initializer, "[[", "strata"))]
            }

            # list for initial compartment volume objects
            for (t in seq_len(n_strata)) {
                comp_prior <-
                    if (!initializer[[t]]$fixed) {
                        if (!is.null(initializer[[t]]$prior)) {
                            if (comp_size_vec[t] != 0) {
                                initializer[[t]]$prior
                            } else {
                                rep(0.0,
                                    length(initializer[[t]]$init_states))
                            }
                        } else {
                            if (comp_size_vec[t] != 0) {
                                initializer[[t]]$init_states
                            } else {
                                rep(0.0,
                                    length(initializer[[t]]$init_states))
                            }
                        }
                    } else {
                        if (comp_size_vec[t] != 0) {
                            initializer[[t]]$init_states
                        } else {
                            rep(0.0, length(initializer[[t]]$init_states))
                        }
                    }

                # compartment probabilities
                comp_probs <- if (sum(comp_prior) != 0) {
                    comp_prior / sum(comp_prior)
                } else {
                    rep(0.0, length(comp_prior))
                }

                # unconstrained moments
                comp_mean  <- comp_size_vec[t] * comp_probs
                comp_cov   <-
                    comp_size_vec[t] * (diag(comp_probs) - comp_probs %*% t(comp_probs))

                if (initializer[[t]]$dist == "dirmultinom")
                    comp_cov <-
                    comp_cov * ((comp_size_vec[t] + sum(comp_prior)) / (1 + sum(comp_prior)))

                comp_cov_svd <- svd(comp_cov)
                comp_cov_svd$d[length(comp_cov_svd$d)] <- 0
                comp_sqrt_cov <-
                    comp_cov_svd$u %*% diag(sqrt(comp_cov_svd$d))

                initdist_objects[[t]] <-
                    list(
                        fixed              = initializer[[t]]$fixed,
                        comp_size          = comp_size_vec[t],
                        comp_mean          = comp_mean,
                        comp_sqrt_cov      = comp_sqrt_cov[, -length(comp_mean)],
                        draws_cur          = rep(0.0, length(comp_mean) - 1),
                        draws_prop         = rep(0.0, length(comp_mean) - 1),
                        draws_ess          = rep(0.0, length(comp_mean) - 1),
                        comp_inds_R        = initializer[[t]]$codes,
                        comp_inds_Cpp      = initializer[[t]]$codes - 1
                    )
            }

            # matrix of initial states
            init_states <-
                matrix(0.0,
                       nrow = nsim,
                       ncol = length(stem_object$dynamics$comp_codes))
            colnames(init_states) <-
                names(stem_object$dynamics$comp_codes)

            for (s in seq_len(stem_object$dynamics$n_strata)) {
                init_states[, stem_object$dynamics$initializer[[s]]$codes] <-
                    matrix(
                        as.numeric(stem_object$dynamics$initializer[[s]]$init_states),
                        nrow = nsim,
                        ncol = length(stem_object$dynamics$initializer[[s]]$init_states),
                        byrow = TRUE
                    )
            }

            # if the initial state is not fixed, sample the collection of initial states
            if (!stem_object$dynamics$fixed_inits) {
                for (n in seq_len(nsim)) {
                    for (s in seq_along(initdist_objects)) {
                        if (!initdist_objects[[s]]$fixed) {
                            # N(0,1) draws
                            draw_normals(initdist_objects[[s]]$draws_cur)

                            # map to volumes
                            copy_vec2(
                                dest = init_state,
                                orig = initdist_objects[[s]]$comp_mean +
                                    c(
                                        initdist_objects[[s]]$comp_sqrt_cov %*% initdist_objects[[s]]$draws_cur
                                    ),
                                inds = initdist_objects[[s]]$comp_inds_Cpp
                            )

                            while (any(init_state[initdist_objects[[s]]$comp_inds_R] < 0) |
                                   any(init_state[initdist_objects[[s]]$comp_inds_R] > initdist_objects[[s]]$comp_size)) {
                                # N(0,1) draws
                                draw_normals(initdist_objects[[s]]$draws_cur)

                                # map to volumes
                                copy_vec2(
                                    dest = init_state,
                                    orig = initdist_objects[[s]]$comp_mean +
                                        c(
                                            initdist_objects[[s]]$comp_sqrt_cov %*% initdist_objects[[s]]$draws_cur
                                        ),
                                    inds = initdist_objects[[s]]$comp_inds_Cpp
                                )
                            }
                        }
                    }

                    # save to initdist matrix
                    init_states[n,] <- c(init_state)
                }

            } else {
                # if all initial states are fixed, just copy the initial compartment counts
                init_states <- matrix(rep(
                    as.numeric(stem_object$dynamics$initdist_params),
                    nsim
                ),
                nrow = nsim,
                byrow = TRUE)
                colnames(init_states) <-
                    names(stem_object$dynamics$comp_codes)
            }

            # tparam indices and initial values
            if (!is.null(stem_object$dynamics$tparam)) {
                # generate indices for time-varying parameters
                for (s in seq_along(stem_object$dynamics$tparam)) {
                    stem_object$dynamics$tparam[[s]]$col_ind  <-
                        stem_object$dynamics$lna_rates$lna_param_codes[stem_object$dynamics$tparam[[s]]$tparam_name]
                    stem_object$dynamics$tparam[[s]]$tpar_inds <-
                        findInterval(lna_times,
                                     stem_object$dynamics$tparam[[s]]$times,
                                     left.open = F) - 1
                    stem_object$dynamics$tparam[[s]]$tpar_inds[stem_object$dynamics$tparam[[s]]$tpar_inds == -1] <-
                        0
                }

                # list for saving the time varying parameters for reuse in simulating a dataset in necessary
                tparam_times <-
                    sort(unique(unlist(
                        lapply(stem_object$dynamics$tparam, function(x)
                            x$times)
                    )))
                tparam_times <-
                    tparam_times[tparam_times >= t0 & tparam_times <= tmax]

                if (is.null(tparam_values) &
                    is.null(tparam_draws)) {
                    tpar_list <-
                        list(lapply(stem_object$dynamics$tparam, function(x)
                            numeric(length(x$values))))

                    tparam_draws  <- rep(tpar_list, nsim)
                    tparam_values <- rep(tpar_list, nsim)

                    # compute the tparam values
                    for (n in seq_len(nsim)) {
                        for (m in seq_along(stem_object$dynamics$tparam)) {
                            # grab parameters
                            if (!is.null(simulation_parameters)) {
                                sim_pars <- simulation_parameters[[n]]
                                class(sim_pars) <- "numeric"
                                lna_pars[lna_param_inds + 1, ] <-
                                    sim_pars
                            }

                            # draw values
                            draw_normals(tparam_draws[[n]][[m]])

                            # map to parameter
                            insert_tparam(
                                tcovar    = lna_pars,
                                values    =
                                    tparam[[m]]$draws2par(
                                        parameters = lna_pars[1, ],
                                        draws = tparam_draws[[n]][[m]]
                                    ),
                                col_ind   = stem_object$dynamics$tparam[[m]]$col_ind,
                                tpar_inds = stem_object$dynamics$tparam[[m]]$tpar_inds
                            )

                            # copy_values
                            tparam_values[[n]][[m]] <-
                                tparam[[m]]$draws2par(parameters = lna_pars[1, ],
                                                      draws = tparam_draws[[n]][[m]])
                        }
                    }

                } else if (is.null(tparam_values) &
                           !is.null(tparam_draws)) {
                    tpar_list <-
                        list(lapply(stem_object$dynamics$tparam, function(x)
                            numeric(length(x$values))))

                    tparam_values <- rep(tpar_list, nsim)

                    # compute the tparam values
                    for (n in seq_len(nsim)) {
                        for (m in seq_along(stem_object$dynamics$tparam)) {
                            if (!is.null(simulation_parameters)) {
                                sim_pars <- simulation_parameters[[n]]
                                class(sim_pars) <- "numeric"
                                lna_pars[lna_param_inds + 1, ] <-
                                    sim_pars
                            }

                            if (!stem_object$dynamics$fixed_inits) {
                                init_vols <- init_states[n, ]
                            }

                            # set the parameters and initial volumes
                            pars2lnapars2(lna_pars, as.numeric(c(
                                sim_pars, init_vols
                            )), 0)

                            # map to parameter
                            insert_tparam(
                                tcovar    = lna_pars,
                                values    =
                                    tparam[[m]]$draws2par(
                                        parameters = lna_pars[1, ],
                                        draws = tparam_draws[[n]][[m]]
                                    ),
                                col_ind   = stem_object$dynamics$tparam[[m]]$col_ind,
                                tpar_inds = stem_object$dynamics$tparam[[m]]$tpar_inds
                            )

                            # copy_values
                            tparam_values[[n]][[m]] <-
                                tparam[[m]]$draws2par(parameters = lna_pars[1, ],
                                                      draws = tparam_draws[[n]][[m]])
                        }
                    }
                }

            } else {
                tparam_draws <- NULL
                tparam_values <- NULL
                tparam_times <- NULL
            }

            # generate forcing matrix
            if (!is.null(stem_object$dynamics$dynamics_args$forcings)) {
                # names and indices
                forcing_tcovars   <-
                    sapply(forcings, function(x)
                        x$tcovar_name)
                forcing_tcov_inds <-
                    match(forcing_tcovars, colnames(lna_pars)) - 1
                forcing_events    <-
                    c(sapply(forcings, function(x)
                        paste0(x$from, "2", x$to)))

                # matrix indicating which compartments are involved in which forcings in and out
                forcings_out <- matrix(
                    0.0,
                    nrow = ncol(stem_object$dynamics$flow_matrix_lna),
                    ncol = length(forcings),
                    dimnames = list(
                        colnames(stem_object$dynamics$flow_matrix_lna),
                        forcing_tcovars
                    )
                )

                forcing_transfers <- array(
                    0.0,
                    dim = c(
                        ncol(stem_object$dynamics$flow_matrix_lna),
                        ncol(stem_object$dynamics$flow_matrix_lna),
                        length(forcings)
                    ),
                    dimnames = list(
                        colnames(stem_object$dynamics$flow_matrix_lna),
                        colnames(stem_object$dynamics$flow_matrix_lna),
                        forcing_tcovars
                    )
                )

                for (s in seq_along(forcings)) {
                    forcings_out[forcings[[s]]$from, s] <- 1

                    for (t in seq_along(forcings[[s]]$from)) {
                        forcing_transfers[forcings[[s]]$from[t], forcings[[s]]$from[t], s] <-
                            -1
                        forcing_transfers[forcings[[s]]$to[t], forcings[[s]]$from[t], s]    <-
                            1
                    }
                }

            } else {
                forcing_tcovars   <- character(0L)
                forcing_tcov_inds <- integer(0L)
                forcing_events    <- character(0L)
                forcings_out      <- matrix(0.0, nrow = 0, ncol = 0)
                forcing_transfers <- array(0.0, dim = c(0, 0, 0))
            }

            # generate some auxilliary objects
            param_update_inds <- lna_times %in%
                sort(unique(
                    c(t0, tmax, stem_object$dynamics$tcovar[, 1], tparam_times)
                ))
            census_interval_inds <-
                findInterval(lna_times, census_times, left.open = T)

            # initialize the volumes and prevalence indices
            init_vols <- init_states[1, ]
            prev_inds <-
                match(round(census_times, digits = 8),
                      round(lna_times, digits = 8))

            for (k in seq_len(nsim)) {
                if (!is.null(simulation_parameters)) {
                    sim_pars <- simulation_parameters[[k]]
                    class(sim_pars) <- "numeric"
                    lna_pars[lna_param_inds + 1, ] <- sim_pars
                }

                if (!stem_object$dynamics$fixed_inits) {
                    init_vols <- init_states[k, ]
                }

                # set the parameters and initial volumes
                pars2lnapars2(lna_pars, as.numeric(c(sim_pars, init_vols)), 0)

                # draw values for the time-varying parameters
                if (!is.null(stem_object$dynamics$tparam)) {
                    for (s in seq_along(stem_object$dynamics$tparam)) {
                        # insert the new values into the tcovar matrix
                        insert_tparam(
                            tcovar    = lna_pars,
                            values    = tparam_values[[k]][[s]],
                            col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                            tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds
                        )
                    }
                }

                attempt <- 0
                path    <- NULL
                while (is.null(path) && (attempt < max_attempts)) {
                    if (lna_method == "exact") {
                        try({
                            path <- propose_lna(
                                lna_times         = lna_census_times,
                                lna_draws         = lna_draws[[k]],
                                lna_pars          = lna_pars,
                                init_start        = stem_object$dynamics$lna_initdist_inds[1],
                                lna_param_inds    = parameter_inds,
                                lna_tcovar_inds   = tcovar_inds,
                                param_update_inds = param_update_inds,
                                stoich_matrix     = stem_object$dynamics$stoich_matrix_lna,
                                forcing_inds      = forcing_inds,
                                forcing_tcov_inds = forcing_tcov_inds,
                                forcings_out      = forcings_out,
                                forcing_transfers = forcing_transfers,
                                step_size         = stem_object$dynamics$dynamics_args$step_size,
                                max_attempts      = max_attempts,
                                lna_pointer       = stem_object$dynamics$lna_pointers$lna_ptr,
                                set_pars_pointer  = stem_object$dynamics$lna_pointers$set_lna_params_ptr
                            )

                        }, silent = TRUE)

                        attempt           <- attempt + 1

                        if (!is.null(path)) {
                            census_paths[[k]] <-
                                census_incidence(path$lna_path,
                                                 census_times,
                                                 census_interval_inds)
                            lna_paths[[k]]    <-
                                path$prev_path[prev_inds, ]
                            lna_draws[[k]]    <- path$draws
                        } else {
                            lna_draws[[k]]    <-
                                matrix(rnorm(lna_draws[[k]]), nrow(lna_draws[[k]]))
                        }

                    } else if (lna_method == "approx") {
                        try({
                            path <- propose_lna_approx(
                                lna_times         = lna_census_times,
                                lna_draws         = lna_draws[[k]],
                                lna_pars          = lna_pars,
                                init_start        = stem_object$dynamics$lna_initdist_inds[1],
                                lna_param_inds    = parameter_inds,
                                lna_tcovar_inds   = tcovar_inds,
                                param_update_inds = param_update_inds,
                                stoich_matrix     = stem_object$dynamics$stoich_matrix_lna,
                                forcing_inds      = forcing_inds,
                                forcing_tcov_inds = forcing_tcov_inds,
                                forcings_out      = forcings_out,
                                forcing_transfers = forcing_transfers,
                                step_size         = stem_object$dynamics$dynamics_args$step_size,
                                max_attempts      = max_attempts,
                                ess_updates       = 1,
                                ess_warmup        = ess_warmup,
                                lna_bracket_width = lna_bracket_width,
                                lna_pointer       = stem_object$dynamics$lna_pointers$lna_ptr,
                                set_pars_pointer  = stem_object$dynamics$lna_pointers$set_lna_params_ptr
                            )
                        }, silent = TRUE)

                        attempt           <- attempt + 1

                        if (!is.null(path)) {
                            census_paths[[k]] <-
                                census_incidence(path$incid_paths,
                                                 census_times,
                                                 census_interval_inds)
                            lna_paths[[k]]    <-
                                path$prev_paths[prev_inds, ]
                            lna_draws[[k]]    <- path$draws
                        } else {
                            lna_draws[[k]]    <-
                                matrix(rnorm(lna_draws[[k]]), nrow(lna_draws[[k]]))
                        }
                    }
                }

                if (is.null(census_paths[[k]])) {
                    warning(
                        "Simulation failed. Increase max attempts per simulation or try different parameter values."
                    )
                } else {
                    colnames(census_paths[[k]]) <-
                        c("time",
                          rownames(stem_object$dynamics$flow_matrix_lna))
                    colnames(lna_paths[[k]]) <-
                        c("time",
                          colnames(stem_object$dynamics$flow_matrix_lna))
                }
            }

            failed_runs  <- which(sapply(lna_paths, is.null))
            if (length(failed_runs) != 0) {
                census_paths <- census_paths[-failed_runs]
                lna_paths    <- lna_paths[-failed_runs]
                lna_draws    <- lna_draws[-failed_runs]
            }

        } else if (method == "ode") {

            # set the vectors of times when the ODE is evaluated and censused
            ode_times <-
                sort(unique(
                    c(t0, census_times, stem_object$dynamics$tcovar[, 1], tmax)
                ))

            ode_census_times <-
                ode_times[ode_times >= t0 & ode_times <= tmax]

            # generate the matrix of parameters, constants, and time-varying covariates
            ode_pars  <- matrix(
                0.0,
                nrow = length(ode_times),
                ncol = length(stem_object$dynamics$ode_rates$ode_param_codes),
                dimnames = list(
                    NULL,
                    names(stem_object$dynamics$ode_rates$ode_param_codes)
                )
            )

            # insert parameters, constants, and time-varying covariates
            parameter_inds <-
                setdiff(
                    stem_object$dynamics$param_codes,
                    stem_object$dynamics$ode_initdist_inds
                )

            constant_inds  <-
                length(stem_object$dynamics$param_codes) + seq_along(stem_object$dynamics$const_codes) - 1

            tcovar_inds    <-
                length(stem_object$dynamics$param_codes) + length(constant_inds) + seq_along(stem_object$dynamics$tcovar_codes) - 1

            ode_pars[, parameter_inds + 1] <- matrix(0.0,
                                                     nrow = nrow(ode_pars),
                                                     ncol = length(parameter_inds))
            ode_pars[, constant_inds + 1]  <- matrix(0.0,
                                                     nrow = nrow(ode_pars),
                                                     ncol = length(constant_inds))

            pars2lnapars2(
                lnapars = ode_pars,
                parameters = c(as.numeric(
                    stem_object$dynamics$parameters[parameter_inds + 1]
                )),
                c_start = 0
            )

            pars2lnapars2(
                lnapars = ode_pars,
                parameters = as.numeric(stem_object$dynamics$constants),
                c_start = constant_inds[1]
            )

            # Generate forcing indices and update the ODE parameter matrix with forcings
            forcing_inds <- rep(FALSE, length(ode_census_times))

            if (!is.null(stem_object$dynamics$dynamics_args$tcovar)) {
                tcovar_rowinds <- match(
                    round(ode_times, digits = 8),
                    round(stem_object$dynamics$tcovar[, 1], digits = 8)
                )
                ode_pars[tcovar_rowinds, tcovar_inds + 1] <-
                    stem_object$dynamics$tcovar[tcovar_rowinds, -1]

                # zero out forcings if necessary
                if (!is.null(stem_object$dynamics$dynamics_args$forcings)) {
                    # get the forcing indices (supplied in the original tcovar matrix)
                    for (f in seq_along(stem_object$dynamics$forcings)) {
                        forcing_inds <- forcing_inds |
                            stem_object$dynamics$tcovar[stem_object$dynamics$tcovar[, 1] %in% ode_census_times,
                                                        stem_object$dynamics$forcings[[f]]$tcovar_name] != 0
                    }

                    zero_inds    <- !forcing_inds

                    # zero out the tcovar elements corresponding to times with no forcings
                    for (l in seq_along(stem_object$dynamics$dynamics_args$forcings)) {
                        ode_pars[zero_inds, stem_object$dynamics$dynamics_args$forcings[[l]]$tcovar_name] = 0
                    }
                }
            }

            # get parameters
            sim_pars  <-
                stem_object$dynamics$parameters[parameter_inds + 1]

            # generate some auxilliary objects
            param_update_inds <-
                ode_times %in% unique(c(t0, tmax, stem_object$dynamics$tcovar[, 1]))

            census_interval_inds <-
                findInterval(ode_times, census_times, left.open = T)

            # objects for storing the ODE paths
            fixed_parameters <-
                stem_object$dynamics$fixed_inits & is.null(simulation_parameters)

            if (!fixed_parameters) {
                census_paths <- vector(mode = "list", length = nsim)
                ode_paths    <- vector(mode = "list", length = nsim)
            } else {
                census_paths <- vector(mode = "list", length = 1)
                ode_paths    <- vector(mode = "list", length = 1)
            }

            # retrieve the initial state
            init_state <-
                ode_pars[1, stem_object$dynamics$ode_initdist_inds + 1]

            # object for sampling initial states
            n_strata <- stem_object$dynamics$n_strata
            initializer <- stem_object$dynamics$initializer
            initdist_objects <-
                vector("list", length = stem_object$dynamics$n_strata)

            if (n_strata == 1) {
                comp_size_vec <- stem_object$dynamics$constants["popsize"]
            } else {
                comp_size_vec <-
                    stem_object$dynamics$constants[paste0("popsize_", sapply(initializer, "[[", "strata"))]
            }

            initdist_objects <-
                vector("list", length = stem_object$dynamics$n_strata)

            if (n_strata == 1) {
                comp_size_vec <- stem_object$dynamics$constants["popsize"]
            } else {
                comp_size_vec <-
                    stem_object$dynamics$constants[paste0("popsize_", sapply(initializer, "[[", "strata"))]
            }

            # list for initial compartment volume objects
            for (t in seq_len(n_strata)) {
                comp_prior <-
                    if (!initializer[[t]]$fixed) {
                        if (!is.null(initializer[[t]]$prior)) {
                            if (comp_size_vec[t] != 0) {
                                initializer[[t]]$prior
                            } else {
                                rep(0.0,
                                    length(initializer[[t]]$init_states))
                            }
                        } else {
                            if (comp_size_vec[t] != 0) {
                                initializer[[t]]$init_states
                            } else {
                                rep(0.0,
                                    length(initializer[[t]]$init_states))
                            }
                        }
                    } else {
                        if (comp_size_vec[t] != 0) {
                            initializer[[t]]$init_states
                        } else {
                            rep(0.0, length(initializer[[t]]$init_states))
                        }
                    }

                # compartment probabilities
                comp_probs <- if (sum(comp_prior) != 0) {
                    comp_prior / sum(comp_prior)
                } else {
                    rep(0.0, length(comp_prior))
                }

                # unconstrained moments
                comp_mean  <- comp_size_vec[t] * comp_probs
                comp_cov   <-
                    comp_size_vec[t] * (diag(comp_probs) - comp_probs %*% t(comp_probs))

                if (initializer[[t]]$dist == "dirmultinom")
                    comp_cov <-
                    comp_cov * ((comp_size_vec[t] + sum(comp_prior)) / (1 + sum(comp_prior)))

                comp_cov_svd <- svd(comp_cov)
                comp_cov_svd$d[length(comp_cov_svd$d)] <- 0
                comp_sqrt_cov <-
                    comp_cov_svd$u %*% diag(sqrt(comp_cov_svd$d))

                initdist_objects[[t]] <-
                    list(
                        fixed              = initializer[[t]]$fixed,
                        comp_size          = comp_size_vec[t],
                        comp_mean          = comp_mean,
                        comp_sqrt_cov      = comp_sqrt_cov[, -length(comp_mean)],
                        draws_cur          = rep(0.0, length(comp_mean) - 1),
                        draws_prop         = rep(0.0, length(comp_mean) - 1),
                        draws_ess          = rep(0.0, length(comp_mean) - 1),
                        comp_inds_R        = initializer[[t]]$codes,
                        comp_inds_Cpp      = initializer[[t]]$codes - 1
                    )
            }

            # matrix of initial states
            init_states <-
                matrix(0.0,
                       nrow = nsim,
                       ncol = length(stem_object$dynamics$comp_codes))

            colnames(init_states) <-
                names(stem_object$dynamics$comp_codes)

            for (s in seq_len(stem_object$dynamics$n_strata)) {
                init_states[, stem_object$dynamics$initializer[[s]]$codes] <-
                    matrix(
                        as.numeric(stem_object$dynamics$initializer[[s]]$init_states),
                        nrow = nsim,
                        ncol = length(stem_object$dynamics$initializer[[s]]$init_states),
                        byrow = TRUE
                    )
            }

            # if the initial state is not fixed, sample the collection of initial states
            if (!stem_object$dynamics$fixed_inits) {
                for (n in seq_len(nsim)) {
                    for (s in seq_along(initdist_objects)) {
                        if (!initdist_objects[[s]]$fixed) {
                            # N(0,1) draws
                            draw_normals(initdist_objects[[s]]$draws_cur)

                            # map to volumes
                            copy_vec2(
                                dest = init_state,
                                orig = initdist_objects[[s]]$comp_mean +
                                    c(
                                        initdist_objects[[s]]$comp_sqrt_cov %*% initdist_objects[[s]]$draws_cur
                                    ),
                                inds = initdist_objects[[s]]$comp_inds_Cpp
                            )

                            while (any(init_state[initdist_objects[[s]]$comp_inds_R] < 0) |
                                   any(init_state[initdist_objects[[s]]$comp_inds_R] > initdist_objects[[s]]$comp_size)) {
                                # N(0,1) draws
                                draw_normals(initdist_objects[[s]]$draws_cur)

                                # map to volumes
                                copy_vec2(
                                    dest = init_state,
                                    orig = initdist_objects[[s]]$comp_mean +
                                        c(
                                            initdist_objects[[s]]$comp_sqrt_cov %*% initdist_objects[[s]]$draws_cur
                                        ),
                                    inds = initdist_objects[[s]]$comp_inds_Cpp
                                )
                            }
                        }
                    }

                    # save to initdist matrix
                    init_states[n,] <- c(init_state)
                }
            }

            # tparam indices and initial values
            if (!is.null(stem_object$dynamics$tparam)) {
                # generate indices for time-varying parameters
                for (s in seq_along(stem_object$dynamics$tparam)) {
                    stem_object$dynamics$tparam[[s]]$col_ind  <-
                        stem_object$dynamics$ode_rates$ode_param_codes[stem_object$dynamics$tparam[[s]]$tparam_name]
                    stem_object$dynamics$tparam[[s]]$tpar_inds <-
                        findInterval(ode_times,
                                     stem_object$dynamics$tparam[[s]]$times,
                                     left.open = F) - 1
                    stem_object$dynamics$tparam[[s]]$tpar_inds[stem_object$dynamics$tparam[[s]]$tpar_inds == -1] <-
                        0
                }

                # list for saving the time varying parameters for reuse in simulating a dataset in necessary
                tparam_times <-
                    sort(unique(unlist(
                        lapply(stem_object$dynamics$tparam, function(x)
                            x$times)
                    )))
                tparam_census_times <-
                    tparam_times[tparam_times >= t0 & tparam_times <= tmax]

                if (is.null(tparam_values) &
                    is.null(tparam_draws)) {
                    tpar_list <-
                        list(lapply(stem_object$dynamics$tparam, function(x)
                            numeric(length(x$values))))

                    tparam_draws  <- rep(tpar_list, nsim)
                    tparam_values <- rep(tpar_list, nsim)

                    # compute the tparam values
                    for (n in seq_len(nsim)) {
                        for (m in seq_along(stem_object$dynamics$tparam)) {
                            if (!is.null(simulation_parameters)) {
                                sim_pars <- simulation_parameters[[n]]
                                class(sim_pars) <- "numeric"
                            }

                            if (!stem_object$dynamics$fixed_inits) {
                                init_vols <- init_states[n, ]
                            }

                            # set the parameters and initial volumes
                            pars2lnapars2(
                                lnapars = ode_pars,
                                parameters = as.numeric(c(
                                    sim_pars, init_vols
                                )),
                                c_start = 0
                            )

                            # draw values
                            draw_normals(tparam_draws[[n]][[m]])

                            # map to parameter
                            insert_tparam(
                                tcovar    = ode_pars,
                                values    =
                                    stem_object$dynamics$tparam[[m]]$draws2par(
                                        parameters = ode_pars[1, ],
                                        draws = tparam_draws[[n]][[m]]
                                    ),
                                col_ind   = stem_object$dynamics$tparam[[m]]$col_ind,
                                tpar_inds = stem_object$dynamics$tparam[[m]]$tpar_inds
                            )

                            # copy_values
                            tparam_values[[n]][[m]] <-
                                stem_object$dynamics$tparam[[m]]$draws2par(
                                    parameters = ode_pars[1, ],
                                    draws = tparam_draws[[n]][[m]])
                        }
                    }

                } else if (is.null(tparam_values) &
                           !is.null(tparam_draws)) {
                    tpar_list <-
                        list(lapply(stem_object$dynamics$tparam, function(x)
                            numeric(length(x$values))))

                    tparam_values <- rep(tpar_list, nsim)

                    # compute the tparam values
                    for (n in seq_len(nsim)) {
                        for (m in seq_along(stem_object$dynamics$tparam)) {
                            # grab parameters
                            if (!is.null(simulation_parameters)) {
                                sim_pars <- simulation_parameters[[n]]
                                class(sim_pars) <- "numeric"
                            }

                            # compute values
                            tparam_values[[n]][[m]] <-
                                stem_object$dynamics$tparam[[m]]$draws2par(sim_pars,
                                                                           tparam_draws[[n]][[m]])
                        }
                    }
                }

            } else {
                tparam_draws <- NULL
                tparam_times <- NULL
            }

            # generate forcing matrix
            if (!is.null(stem_object$dynamics$dynamics_args$forcings)) {
                # names and indices
                forcing_tcovars   <-
                    sapply(forcings, function(x)
                        x$tcovar_name)
                forcing_tcov_inds <-
                    match(forcing_tcovars, colnames(ode_pars)) - 1
                forcing_events    <-
                    c(sapply(forcings, function(x)
                        paste0(x$from, "2", x$to)))

                # matrix indicating which compartments are involved in which forcings in and out
                forcings_out <- matrix(
                    0.0,
                    nrow = ncol(stem_object$dynamics$flow_matrix_ode),
                    ncol = length(forcings),
                    dimnames = list(
                        colnames(stem_object$dynamics$flow_matrix_ode),
                        forcing_tcovars
                    )
                )

                forcing_transfers <- array(
                    0.0,
                    dim = c(
                        ncol(stem_object$dynamics$flow_matrix_ode),
                        ncol(stem_object$dynamics$flow_matrix_ode),
                        length(forcings)
                    ),
                    dimnames = list(
                        colnames(stem_object$dynamics$flow_matrix_ode),
                        colnames(stem_object$dynamics$flow_matrix_ode),
                        forcing_tcovars
                    )
                )

                for (s in seq_along(forcings)) {
                    forcings_out[forcings[[s]]$from, s] <- 1

                    for (t in seq_along(forcings[[s]]$from)) {
                        forcing_transfers[forcings[[s]]$from[t], forcings[[s]]$from[t], s] <-
                            -1
                        forcing_transfers[forcings[[s]]$to[t], forcings[[s]]$from[t], s]    <-
                            1
                    }
                }

            } else {
                forcing_tcovars   <- character(0L)
                forcing_tcov_inds <- integer(0L)
                forcing_events    <- character(0L)
                forcings_out      <- matrix(0.0, nrow = 0, ncol = 0)
                forcing_transfers <- array(0.0, dim = c(0, 0, 0))
            }

            sim_pars  <-
                stem_object$dynamics$parameters[parameter_inds + 1]
            init_vols <- init_states[1, ]
            prev_inds <-
                match(round(census_times, digits = 8),
                      round(ode_times, digits = 8))

            for (k in seq_along(census_paths)) {
                if (!is.null(simulation_parameters)) {
                    sim_pars <- simulation_parameters[[k]]
                    class(sim_pars) <- "numeric"
                }

                if (!stem_object$dynamics$fixed_inits) {
                    init_vols <- init_states[k, ]
                }

                # set the parameters and initial volumes
                pars2lnapars2(
                    lnapars = ode_pars,
                    parameters = as.numeric(c(sim_pars, init_vols)),
                    c_start = 0
                )

                # draw values for the time-varying parameters
                if (!is.null(stem_object$dynamics$tparam)) {
                    for (s in seq_along(stem_object$dynamics$tparam)) {
                        # insert the new values into the tcovar matrix
                        insert_tparam(
                            tcovar    = ode_pars,
                            values    = tparam_values[[k]][[s]],
                            col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                            tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds
                        )
                    }
                }

                path <- NULL

                try({
                    path <- integrate_odes(
                        ode_times         = ode_times,
                        ode_pars          = ode_pars,
                        init_start        = stem_object$dynamics$ode_initdist_inds[1],
                        ode_param_inds    = parameter_inds,
                        ode_tcovar_inds   = tcovar_inds,
                        param_update_inds = param_update_inds,
                        stoich_matrix     = stem_object$dynamics$stoich_matrix_ode,
                        forcing_inds      = forcing_inds,
                        forcing_tcov_inds = forcing_tcov_inds,
                        forcings_out      = forcings_out,
                        forcing_transfers = forcing_transfers,
                        step_size         = stem_object$dynamics$dynamics_args$step_size,
                        ode_pointer       = stem_object$dynamics$ode_pointers$ode_ptr,
                        set_pars_pointer  = stem_object$dynamics$ode_pointers$set_ode_params_ptr
                    )
                }, silent = TRUE)

                if (!is.null(path)) {
                    census_paths[[k]] <-
                        census_incidence(path$incid_path,
                                         census_times,
                                         census_interval_inds)
                    ode_paths[[k]]    <-
                        path$prev_path[match(round(census_times, digits = 8),
                                             round(path$prev_path[, 1], digits = 8)), ]

                    colnames(census_paths[[k]]) <-
                        c("time",
                          rownames(stem_object$dynamics$flow_matrix_ode))
                    colnames(ode_paths[[k]]) <-
                        c("time",
                          colnames(stem_object$dynamics$flow_matrix_ode))
                }

                if (is.null(path) & messages) {
                    warning("Simulation failed. Try different parameter values.")
                }
            }

            failed_runs  <- which(sapply(ode_paths, is.null))
            if (length(failed_runs) != 0) {
                census_paths <- census_paths[-failed_runs]
                ode_paths    <- ode_paths[-failed_runs]
            }
        }

        if (observations && length(failed_runs) != nsim) {

            datasets <- vector(mode = "list", length = length(census_paths))
            measvar_names <- colnames(stem_object$measurement_process$obsmat)

            # grab the time-varying covariate values at observation times
            tcovar_obstimes <-
                build_census_path(
                    path           = stem_object$dynamics$tcovar,
                    census_times   = stem_object$measurement_process$obstimes,
                    census_columns = 1:(ncol(stem_object$dynamics$tcovar) - 1)
                )

            colnames(tcovar_obstimes) <-
                colnames(stem_object$dynamics$tcovar)

            if (method != "gillespie") {
                # ditch the first column of tcovar_obstimes
                tcovar_obstimes = tcovar_obstimes[, -1, drop = FALSE]
            }

            # if incidence, the incidence codes are not null
            do_incidence <-
                !is.null(stem_object$dynamics$incidence_codes) &
                !(method %in% c("lna", "ode"))

            # census if computing incidence or if computing prevalence and the obstimes != census_times
            do_census <- !is.null(census_times) &&
                !identical(
                    as.numeric(stem_object$measurement_process$obstimes),
                    as.numeric(census_times)
                )

            # get the indices in the censused matrices for the observation times
            if (do_census)
                cens_inds <-
                findInterval(stem_object$measurement_process$obstimes,
                             census_times)

            if (method == "gillespie") {
                measproc_indmat = as.matrix(stem_object$measurement_process$measproc_indmat)
                constants       = as.numeric(stem_object$dynamics$constants)
                r_measure_ptr   = stem_object$measurement_process$meas_pointers$r_measure_ptr

                # should incidence be computed?
                if (do_incidence)
                    incidence_codes <- stem_object$dynamics$incidence_codes + 1

                # column codes in the path matrix for compartments to be censused
                census_codes <-
                    c(stem_object$dynamics$comp_codes,
                      stem_object$dynamics$incidence_codes) + 2

                pathmat         <-
                    stem_object$measurement_process$censusmat

                # initialize simulation parameters
                sim_pars <- stem_object$dynamics$parameters
                class(sim_pars) <- "numeric"

                for (k in seq_len(nsim)) {
                    # get the state at observation times
                    if (!is.null(census_paths[[k]])) {
                        if (do_census) {
                            pathmat[, -1] <- census_paths[[k]][cens_inds,-1]
                        } else {
                            pathmat <- census_paths[[k]]
                        }

                        # get the new simulation parameters if a list was supplied
                        if (!is.null(simulation_parameters)) {
                          sim_pars <- simulation_parameters[[k]]
                          class(sim_pars) <- "numeric"
                        }

                        # insert the time-varying parameters into the tcovar matrix
                        if (!is.null(tparam_draws)) {
                            for (s in seq_along(stem_object$dynamics$tparam)) {
                                # insert the tparam values into the tcovar matrix
                                insert_tparam(
                                    tcovar = stem_object$dynamics$tcovar,
                                    values = tparam_values[[k]][[s]],
                                    col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                                    tpar_inds =
                                        stem_object$dynamics$tparam[[s]]$tpar_inds
                                )
                            }

                            # grab the time-varying covariate values at observation times
                            tcovar_obstimes <-
                                build_census_path(
                                    path           = stem_object$dynamics$tcovar,
                                    census_times = stem_object$measurement_process$obstimes,
                                    census_columns = 1:(
                                        ncol(stem_object$dynamics$tcovar) - 1
                                    )
                                )

                            colnames(tcovar_obstimes) <-
                                colnames(stem_object$dynamics$tcovar)
                        }

                        # simulate the data
                        datasets[[k]] <-
                            simulate_r_measure(
                                censusmat = pathmat,
                                measproc_indmat = measproc_indmat,
                                parameters = sim_pars,
                                constants = constants,
                                tcovar = tcovar_obstimes,
                                r_measure_ptr = r_measure_ptr
                            )
                        colnames(datasets[[k]]) <- measvar_names
                    }
                }

            } else if (method == "lna") {
                # get the objects for simulating from the measurement process
                measproc_indmat  <-
                    stem_object$measurement_process$measproc_indmat
                sim_pars <- stem_object$dynamics$parameters
                class(sim_pars) <- "numeric"
                constants        <-
                    as.numeric(stem_object$dynamics$constants)
                tcovar           <- stem_object$dynamics$tcovar
                r_measure_ptr    <-
                    stem_object$measurement_process$meas_pointers_lna$r_measure_ptr
                cens_inds        <-
                    c(0, match(
                        round(
                            stem_object$measurement_process$obstimes,
                            digits = 8
                        ),
                        round(census_times, digits = 8)
                    ) - 1)
                do_prevalence    <-
                    stem_object$measurement_process$lna_prevalence
                do_incidence     <-
                    stem_object$measurement_process$lna_incidence
                obstime_inds     <-
                    stem_object$measurement_process$obstime_inds
                pathmat          <-
                    stem_object$measurement_process$censusmat
                flow_matrix_lna  <-
                    stem_object$dynamics$flow_matrix_lna
                event_inds   <-
                    stem_object$measurement_process$incidence_codes_lna

                # reinitialize the tparam indices if necessary
                if (!is.null(stem_object$dynamics$tparam)) {
                    # reinitialize the tparam indices if necessary
                    for (s in seq_along(stem_object$dynamics$tparam)) {
                        stem_object$dynamics$tparam[[s]]$col_ind <-
                            stem_object$dynamics$tcovar_codes[stem_object$dynamics$tparam[[s]]$tparam_name]

                        stem_object$dynamics$tparam[[s]]$tpar_inds <-
                            findInterval(
                                stem_object$dynamics$tcovar[, 1],
                                stem_object$dynamics$tparam[[s]]$times,
                                left.open = F
                            ) - 1

                        stem_object$dynamics$tparam[[s]]$tpar_inds[stem_object$dynamics$tparam[[s]]$tpar_inds == -1] <-
                            0
                    }
                }

                # ditch the time column of tcovar
                tcovar <- tcovar[, -1, drop = FALSE]

                for (k in seq_along(census_paths)) {

                    # fill out the census matrix
                    sim_path =
                        cbind(lna_paths[[k]][cens_inds,],
                              census_paths[[k]][cens_inds,-1])[,colnames(pathmat)]

                    # census_latent_path(
                    #     path                = census_paths[[k]],
                    #     census_path         = pathmat,
                    #     census_inds         = cens_inds,
                    #     event_inds          = event_inds,
                    #     flow_matrix         = flow_matrix_lna,
                    #     do_prevalence       = do_prevalence,
                    #     parmat              = ode_pars,
                    #     initdist_inds       = init_states[k, ],
                    #     forcing_inds        = forcing_inds,
                    #     forcing_tcov_inds   = forcing_tcov_inds,
                    #     forcings_out        = forcings_out,
                    #     forcing_transfers   = forcing_transfers
                    # )

                    if (!is.null(simulation_parameters)) {
                      sim_pars <- simulation_parameters[[k]]
                      class(sim_pars) <- "numeric"
                    }


                    # insert the time-varying parameters into the tcovar matrix
                    if (!is.null(tparam_draws)) {
                        for (s in seq_along(stem_object$dynamics$tparam)) {
                            # insert the tparam values into the tcovar matrix
                            insert_tparam(
                                tcovar = stem_object$dynamics$tcovar,
                                values = tparam_values[[k]][[s]],
                                col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                                tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds
                            )
                        }

                        # grab the time-varying covariate values at observation times
                        tcovar_obstimes <-
                            build_census_path(
                                path = stem_object$dynamics$tcovar,
                                census_times = stem_object$measurement_process$obstimes,
                                census_columns = 1:(ncol(
                                    stem_object$dynamics$tcovar
                                ) - 1)
                            )

                        colnames(tcovar_obstimes) <-
                            colnames(stem_object$dynamics$tcovar)

                        # ditch the first column of tcovar_obstimes
                        tcovar_obstimes = tcovar_obstimes[, -1, drop = FALSE]
                    }

                    # simulate the dataset
                    datasets[[k]] <- simulate_r_measure(
                        sim_path,
                        # pathmat,
                        measproc_indmat,
                        sim_pars,
                        constants,
                        tcovar_obstimes,
                        r_measure_ptr
                    )
                    colnames(datasets[[k]]) <- measvar_names
                }

            } else if (method == "ode") {

                # get the objects for simulating from the measurement process
                measproc_indmat  <-
                    stem_object$measurement_process$measproc_indmat
                sim_pars <- stem_object$dynamics$parameters
                class(sim_pars) <- "numeric"
                constants        <-
                    as.numeric(stem_object$dynamics$constants)
                tcovar           <- stem_object$dynamics$tcovar
                r_measure_ptr    <-
                    stem_object$measurement_process$meas_pointers_lna$r_measure_ptr
                cens_inds        <-
                    match(
                        round(stem_object$measurement_process$obstimes, digits = 8),
                        round(census_times, digits = 8))
                do_prevalence    <-
                    stem_object$measurement_process$ode_prevalence
                do_incidence     <-
                    stem_object$measurement_process$ode_incidence
                obstime_inds     <-
                    stem_object$measurement_process$obstime_inds
                pathmat          <-
                    stem_object$measurement_process$censusmat
                flow_matrix_ode  <-
                    stem_object$dynamics$flow_matrix_ode
                event_inds   <-
                    stem_object$measurement_process$incidence_codes_ode

                if (!is.null(stem_object$dynamics$tparam)) {

                    # reinitialize the tparam indices if necessary
                    for (s in seq_along(stem_object$dynamics$tparam)) {
                        stem_object$dynamics$tparam[[s]]$col_ind  <-
                            stem_object$dynamics$tcovar_codes[stem_object$dynamics$tparam[[s]]$tparam_name]
                        stem_object$dynamics$tparam[[s]]$tpar_inds <-
                            findInterval(
                                stem_object$dynamics$tcovar[, 1],
                                stem_object$dynamics$tparam[[s]]$times,
                                left.open = F
                            ) - 1

                        stem_object$dynamics$tparam[[s]]$tpar_inds[stem_object$dynamics$tparam[[s]]$tpar_inds == -1] <-
                            0
                    }
                }

                # ditch the time column of tcovar
                tcovar <- tcovar[, -1, drop = FALSE]

                if (!fixed_parameters) {
                    for (k in seq_along(census_paths)) {

                        sim_path =
                            cbind(ode_paths[[k]][cens_inds,],
                                  census_paths[[k]][cens_inds,-1])[,colnames(pathmat)]

                        # census_latent_path(
                        #     path                = sim_path,
                        #     census_path         = pathmat,
                        #     census_inds         = cens_inds,
                        #     event_inds          = event_inds,
                        #     flow_matrix         = flow_matrix_ode,
                        #     do_prevalence       = do_prevalence,
                        #     parmat              = ode_pars,
                        #     initdist_inds       = init_states[k, ],
                        #     forcing_inds        = forcing_inds,
                        #     forcing_tcov_inds   = forcing_tcov_inds,
                        #     forcings_out        = forcings_out,
                        #     forcing_transfers   = forcing_transfers
                        # )

                        if (!is.null(simulation_parameters)) {
                          sim_pars <- simulation_parameters[[k]]
                          class(sim_pars) <- "numeric"
                        }

                        # insert the time-varying parameters into the tcovar matrix
                        if (!is.null(tparam_values)) {
                            for (s in seq_along(stem_object$dynamics$tparam)) {
                                # insert the tparam values into the tcovar matrix
                                insert_tparam(
                                    tcovar = stem_object$dynamics$tcovar,
                                    values = tparam_values[[k]][[s]],
                                    col_ind =
                                        stem_object$dynamics$tparam[[s]]$col_ind,
                                    tpar_inds =
                                        stem_object$dynamics$tparam[[s]]$tpar_inds
                                )
                            }

                            # grab the time-varying covariate values at observation times
                            tcovar_obstimes <-
                                build_census_path(
                                    path = stem_object$dynamics$tcovar,
                                    census_times =
                                        stem_object$measurement_process$obstimes,

                                    census_columns = 1:(
                                        ncol(stem_object$dynamics$tcovar) - 1
                                    )
                                )
                            colnames(tcovar_obstimes) <-
                                colnames(stem_object$dynamics$tcovar)

                            # ditch the first column of tcovar_obstimes
                            tcovar_obstimes = tcovar_obstimes[, -1, drop = FALSE]
                        }

                        # simulate the dataset
                        datasets[[k]] <- simulate_r_measure(
                            sim_path,
                            # pathmat,
                            measproc_indmat,
                            sim_pars,
                            constants,
                            tcovar_obstimes,
                            r_measure_ptr
                        )
                        colnames(datasets[[k]]) <- measvar_names
                    }
                } else {

                    sim_path =
                        cbind(ode_paths[[1]][cens_inds,],
                              census_paths[[1]][cens_inds,-1])[,colnames(pathmat)]

                    # census_latent_path(
                    #     path                = census_paths[[1]],
                    #     census_path         = pathmat,
                    #     census_inds         = cens_inds,
                    #     event_inds          = event_inds,
                    #     flow_matrix         = flow_matrix_ode,
                    #     do_prevalence       = do_prevalence,
                    #     parmat              = ode_pars,
                    #     initdist_inds       = init_states[1, ],
                    #     forcing_inds        = forcing_inds,
                    #     forcing_tcov_inds   = forcing_tcov_inds,
                    #     forcings_out        = forcings_out,
                    #     forcing_transfers   = forcing_transfers
                    # )

                    for (k in seq_len(nsim)) {
                        # insert the time-varying parameters into the tcovar matrix
                        if (!is.null(tparam_draws)) {
                            for (s in seq_along(stem_object$dynamics$tparam)) {
                                # insert the tparam values into the tcovar matrix
                                insert_tparam(
                                    tcovar = stem_object$dynamics$tcovar,
                                    values = tparam_values[[k]][[s]],
                                    col_ind   = stem_object$dynamics$tparam[[s]]$col_ind,
                                    tpar_inds = stem_object$dynamics$tparam[[s]]$tpar_inds
                                )
                            }

                            # grab the time-varying covariate values at observation times
                            tcovar_obstimes <-
                                build_census_path(
                                    path = stem_object$dynamics$tcovar,
                                    census_times   =
                                        stem_object$measurement_process$obstimes,
                                    census_columns = 1:(
                                        ncol(stem_object$dynamics$tcovar) - 1
                                    )
                                )
                            colnames(tcovar_obstimes) <-
                                colnames(stem_object$dynamics$tcovar)

                            # ditch the first column of tcovar_obstimes
                            tcovar_obstimes = tcovar_obstimes[, -1, drop = FALSE]
                        }

                        # simulate the dataset
                        datasets[[k]] <- simulate_r_measure(
                            sim_path,
                            # pathmat,
                            measproc_indmat,
                            sim_pars,
                            constants,
                            tcovar_obstimes,
                            r_measure_ptr
                        )
                        colnames(datasets[[k]]) <- measvar_names
                    }
                }
            }
        } else if (observations && length(failed_runs == nsim)) {
            datasets = list(NULL)
        }

        stem_simulations <- list(paths = NULL, datasets = NULL)

        if (full_paths)
            stem_simulations$full_paths <- paths_full


        if (paths & !is.null(census_times)) {
            stem_simulations$paths <- census_paths
            if (method == "lna")
                stem_simulations$natural_paths <- lna_paths
            if (method == "ode")
                stem_simulations$natural_paths <- ode_paths

        } else {
            stem_simulations$paths <- NULL
        }

        if (observations)
            stem_simulations$datasets      <- datasets
        if (method == "lna")
            stem_simulations$lna_draws   <- lna_draws
        stem_simulations$failed_runs <- failed_runs

        class(stem_simulations) <- "stemr_simulation_list"
        return(stem_simulations)
    }
