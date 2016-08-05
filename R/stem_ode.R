#' Obtain the deterministic ODE path for a stochastic epidemic model.
#'
#' @param stem_object stochastic epidemic model object
#' @param census_times list of times at which to evaluate the path
#' @param init_states initial state
#'
#' @return matrix with the ODE path
#' @export
stem_ode <- function(stem_object, census_times, init_states) {

        # initialize the path matrix
        ode_path        <- matrix(0, nrow = length(census_times), ncol = 1 + length(init_states)) # create matrix
        ode_path[1,]    <- c(census_times[1], init_states)        # assign initial state values
        ode_path[,1]    <- census_times                           # insert the census times

        # pull the incidence codes and compartment codes
        comp_codes <- stem_object$dynamics$comp_codes + 1

        if(is.null(stem_object$dynamics$incidence_codes)) {
                incidence_codes <- NULL

        } else {
                incidence_codes   <- stem_object$dynamics$incidence_codes + 1
                incidence_sources <- stem_object$dynamics$incidence_sources + 1
        }

        # construct the time-varying covariate matrix
        if(is.null(stem_object$dynamics$.dynamics_args$tcovar)) {
                        ode_tcovar              <- matrix(census_times, nrow = length(census_times), ncol = 2)
                        colnames(ode_tcovar)    <- c("_time", "TIME")

        } else {
                        ode_tcovar_timeseq      <- sort(unique(c(census_times, stem_object$dynamics$.dynamics_args$tcovar[,1])))
                        tcovar_inds             <- findInterval(ode_tcovar_timeseq, stem_object$dynamics$.dynamics_args$tcovar[,1])
                        ode_tcovar              <- matrix(c(ode_tcovar_timeseq, stem_object$dynamics$.dynamics_args$tcovar[tcovar_inds,-1]))
                        colnames(ode_tcovar)    <- c("_time", names(stem_object$dynamics$tcovar_codes))
        }

        # initialize the deterministic process
        ode_state <- init_states

        # initialize the model list
        stemr_odemod    <- list(ode_fcn     = NA,
                                parameters  = stem_object$dynamics$parameters,
                                constants   = stem_object$dynamics$constants,
                                tcovar      = ode_tcovar,
                                flow_matrix = t(stem_object$dynamics$flow_matrix),
                                hazard_ptr  = stem_object$dynamics$lna_ptrs$hazard_ptr)

        # create the ODE function
        ode_fcn <- function(t, ode_state, stemr_odemod = stemr_odemod) {
                return(list(stemr_odemod$flow_matrix %*% COMPUTE_HAZARD(t, ode_state, stemr_odemod$parameters, stemr_odemod$constants, stemr_odemod$tcovar, stemr_odemod$hazard_ptr)))
        }

        # insert the ode function into the odemod list
        stemr_odemod$ode_fcn <- ode_fcn

        for(k in 2:(length(census_times))) {

                # increment the interval endpoints
                t_L <- census_times[k-1]
                t_R <- census_times[k]

                # retrieve the time-varying covariate values
                stemr_odemod$tcovar <- ode_tcovar[k-1, -1, drop = FALSE]

                # solve the odes
                ode_state <- stem_ode_path(stemr_odemod, ode_state, t_L, t_R)

                # insert the new path into the path matrix
                ode_path[k,-1] <- ode_state
        }

        colnames(ode_path) <- c("time", names(stem_object$dynamics$comp_codes), names(stem_object$dynamics$incidence_codes))

        # if some of the compartments track incidence, compute the incidence
        if(!is.null(stem_object$dynamics$incidence_codes)) {
                census_incidence_rows <- rep(list(seq_along(census_times) - 1), length(stem_object$dynamics$incidence_codes))
                compute_incidence(censusmat = ode_path, col_inds  = incidence_codes, row_inds  = census_incidence_rows)
        }

        return(ode_path)
}