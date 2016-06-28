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
#'   or "LNA" if simulating paths via the Linear Noise Approximation.
#' @param t0 the time at which the system is initialized. If not supplied, it is
#'   taken to be the first observation time if that is available, or otherwise
#'   0. If supplied, it must be no later than the first observation time.
#' @param tmax the time at which simulation of the system is terminated. If not
#'   supplied, defaults to the last observation time if not supplied.
#' @param timestep time discretization for smooth time-varying functions,
#'   specifically seasonal terms, to be used in simulation when
#'   method="gillespie". If not supplied, \code{timestep} will default to
#'   \code{(tmax - t0)/50}.
#' @param census_interval the time interval at which compartment counts should
#'   be evaluated. Required for \code{method = "LNA"}. If supplied and
#'   \code{method = "gillespie"}, the compartment counts at census intervals are
#'   returned rather than the full paths.
#' @param as_array if TRUE, the simulated paths and/or simulated datasets will
#'   be returned as an array, or as a list containing multiple arrays (one for
#'   the paths, another for the datasets, and a third for subject-level paths if
#'   they are requested), rather than as a list of paths. Available only if
#'   \code{census_interval = TRUE}.
#'
#' @return If \code{paths = FALSE} and \code{observations = FALSE}, or if
#'   \code{paths = TRUE} and \code{observations = TRUE}, a list or array of
#'   \code{nsim stem} paths and datasets is returned. If \code{subject_paths =
#'   TRUE} and \code{method = "gillespie"}, a list of subject-level paths is
#'   returned.
#'
#'   If \code{paths = TRUE} and \code{observations = FALSE}, a list or array of
#'   simulated population-level paths is returned. If \code{method =
#'   "gillespie"} and \code{subject_paths = TRUE}, a list of two arrays is
#'   returned, with one array for population-level paths, and another for the
#'   subject-level mappings.
#'
#'   If \code{paths = FALSE} and \code{observations = TRUE}, a list or array of
#'   simulated datasets is returned.
#' @export
simulate_stem <- function(stem_object, nsim = 1, paths = FALSE, observations = FALSE, subject_paths = FALSE, method = "gillespie", t0 = NULL, tmax = NULL, timestep = NULL,  census_interval = NULL, as_array = FALSE) {

        # ensure that the method is correctly specified
        if(!method %in% c("gillespie", "LNA")) {
                stop("The simulation method must either be 'gillespie' or 'LNA'.")
        }

        # if LNA, subject paths are not available
        if(!subject_paths && (method == "LNA")) {
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

        # specify the time interval within which the process should be simulated
        if(is.null(t0)) {
                if(is.null(obstimes)) {
                        t0 <- 0
                } else {
                        t0 <- min(sapply(stem_object$measurement_process$obstimes, min))
                }
        }

        if(is.null(tmax)) {
                if(is.null(stem_object$measurement_process$obstimes)) {
                        stop("A stopping time for the simulation must be specified either in the tmax argument, or through the observation times.")
                } else {
                        tmax <- max(sapply(stem_object$measurement_process$obstimes))
                }
        }

        # compute the time discretization interval if it is not supplied and there are smooth time-varying covariates
        if(method == "gillespie" && is.null(timestep) && stem_object$dynamics$timevarying) {
                timestep <- (tmax - t0)/50

        } else if(is.null(timestep) && !stem_object$dynamics$timevarying){
                timestep <- tmax - t0
        }

        # build the time varying covariate matrix (includes, at a minimum, the endpoints of the simulation interval)
        # if timestep is null, there are no time-varying covariates
        if(method == "gillespie") {

                # rebuild the time-varying covariate matrix so that it contains the census intervals
                stem_object$dynamics$tcovar <- build_tcovar_matrix(tcovar = stem_object$dynamics$tcovar, timestep = timestep, t0 = t0, tmax = tmax)
                stem_object$dynamics$tcovar_codes <- seq_len(ncol(stem_object$dynamics$tcovar) - 1)
                names(stem_object$dynamics$tcovar_codes) <- colnames(stem_object$dynamics$tcovar)[2:ncol(stem_object$dynamics$tcovar)]
                stem_object$dynamics$n_tcovar <- ncol(stem_object$dynamics$tcovar) - 1
                stem_object$dynamics$tcovar_changemat <- build_tcovar_changemat(stem_object$dynamics$tcovar)

                timecode <- which(names(stem_object$dynamics$tcovar_codes) == "TIME")

                # make the substitution in the rate functions for TIME
                for(s in seq_along(stem_object$dynamics$rates)) {
                        stem_object$dynamics$rates[[s]]$lumped <- gsub("TIME", paste0("tcovar[",timecode,"]"), stem_object$dynamics$rates[[s]]$lumped)
                        stem_object$dynamics$rates[[s]]$unlumped <- gsub("TIME", paste0("tcovar[",timecode,"]"), stem_object$dynamics$rates[[s]]$unlumped)
                }

                # rebuild the time-varying covariate adjacency matrix
                stem_object$dynamics$tcovar_adjmat <- build_tcovar_adjmat(stem_object$dynamics$rates, stem_object$dynamics$tcovar_codes)

                rate_ptrs <- parse_rates(rates = stem_object$dynamics$rates,
                                param_codes = stem_object$dynamics$param_codes,
                                compartment_codes = stem_object$dynamics$comp_codes,
                                const_codes = stem_object$dynamics$const_codes,
                                tcovar_codes = stem_object$dynamics$tcovar_codes)

                # generate or copy the initial states
                if(stem_object$dynamics$fixed_inits) {
                        # if all initial states are fixed, just copy the initial compartment counts
                        init_states <- matrix(stem_object$dynamics$initdist_params, nrow = 1, byrow = TRUE)
                        colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                } else {
                        if(stem_object$dynamics$n_strata == 1) {
                                init_states <- as.matrix(t(rmultinom(nsim, stem_object$dynamics$popsize, stem_object$dynamics$initdist_params)))
                                colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                        } else if(stem_object$dynamics$n_strata > 1) {
                                # generate the matrix of initial compartment counts
                                init_states <- matrix(0, nrow = nsim, ncol = stem_object$dynamics$n_compartments)
                                colnames(init_states) <- names(stem_object$dynamics$comp_codes)

                                for(s in seq_len(stem_object$dynamics$n_strata)) {
                                        init_states[,stem_object$dynamics$state_initializer[[s]]$param_inds] <- as.matrix(t(rmultinom(nsim, stem_object$dynamics$strata_sizes[s], stem_object$dynamics$state_initializer[[s]]$init_states)))
                                }
                        }
                }

                # initialize the list of paths
                paths <- vector(mode = "list", length = nsim)

                # guess the initial dimensions. need an extra column for event times and another for event IDs.
                if(stem_object$dynamics$progressive & any(stem_object$dynamics$absorbing_states)) {
                        if(stem_object$dynamics$n_strata == 1) {
                                init_dims <- c(n_rows = sum(stem_object$dynamics$popsize * stem_object$dynamics$n_compartments), n_cols = stem_object$dynamics$n_compartments + 2)
                        } else {
                                init_dims <- c(n_rows = sum(stem_object$dynamics$strata_sizes * sapply(sapply(stem_object$dynamics$state_initializer, "[[", 4), length)), n_cols = stem_object$dynamics$n_compartments + 2)

                        }
                } else {
                        if(stem_object$dynamics$n_strata == 1) {
                                init_dims <- c(n_rows = sum(stem_object$dynamics$popsize * stem_object$dynamics$n_compartments) * 3, n_cols = stem_object$dynamics$n_compartments + 2)
                        } else if(stem_object$dynamics$n_strata > 1) {
                                init_dims <- c(n_rows = sum(stem_object$dynamics$strata_sizes * sapply(sapply(stem_object$dynamics$state_initializer, "[[", 4), length)) * 3, n_cols = stem_object$dynamics$n_compartments + 2)
                        }
                }

                # make the initial dimensions a little bigger (round up to nearest power of 2)
                p <- 1
                while(2^p < init_dims[1]) {
                        p <- p+1
                        if(2^p > init_dims[1]) init_dims[1] <- 2^p
                }

                # get the compartment names
                path_colnames <- c("time", "event", names(stem_object$dynamics$comp_codes))

                # simulate the paths
                for(k in seq_len(nsim)) {
                        paths[[k]] <- simulate_gillespie(flow             = stem_object$dynamics$flow_matrix,
                                                         parameters       = stem_object$dynamics$parameters,
                                                         constants        = stem_object$dynamics$constants,
                                                         tcovar           = stem_object$dynamics$tcovar,
                                                         init_states      = init_states[k,],
                                                         rate_adjmat      = stem_object$dynamics$rate_adjmat,
                                                         tcovar_adjmat    = stem_object$dynamics$tcovar_adjmat,
                                                         tcovar_changemat = stem_object$dynamics$tcovar_changemat,
                                                         init_dims        = init_dims,
                                                         rate_ptr         = rate_ptrs[[1]])
                        colnames(paths[[k]]) <- path_colnames
                }


        } else {

        }
}