#' Return the ODE path for a stochastic epidemic model.
#'
#' @param stem_object stochastic epidemic model object list
#' @param init_states vector of initial compartment counts
#' @param census_times times at which to evaluate the ode path
#'
#' @return matrix with the deterministic path evaluated at census times
#' @export
stem_ode <- function(stem_object, init_states, census_times) {

        # set the vectors of times when the LNA is evaluated, restarted, and censused
        ode_times       <- sort(unique(c(census_times, stem_object$dynamics$.dynamics_args$tcovar[,1])))
        census_inds     <- !is.na(match(ode_times, census_times))

        # generate the matrix of parameters, constants, and time-varying covariates
        ode_pars <- matrix(0.0, nrow = length(ode_times), ncol = length(stem_object$dynamics$ode_param_codes))
        colnames(ode_pars) <- c(names(stem_object$dynamics$ode_param_codes))

        # insert parameters, constants, and time-varying covariates
        parameter_inds <- 1:length(stem_object$dynamics$param_codes)
        constant_inds  <- (length(parameter_inds)+1):(length(parameter_inds) + length(stem_object$dynamics$const_codes))
        tcovar_inds  <- (max(constant_inds)+1):ncol(ode_pars)
        ode_pars[, parameter_inds] <- matrix(stem_object$dynamics$parameters, nrow = nrow(ode_pars), ncol = length(parameter_inds), byrow = T)
        ode_pars[, constant_inds]  <- matrix(stem_object$dynamics$constants, nrow = nrow(ode_pars), ncol = length(constant_inds), byrow = T)

        if(!is.null(stem_object$dynamics$.dynamics_args$tcovar)) {
                tcovar_rowinds          <- findInterval(ode_times, stem_object$dynamics$.dynamics_args$tcovar[,1])
                ode_pars[, tcovar_inds] <- stem_object$dynamics$.dynamics_args$tcovar[tcovar_rowinds,-1]
        }

        # if there are artificial incidence compartments, copy the
        # incidence counts and add them to the initial state matrix
        if(!is.null(stem_object$dynamics$incidence_codes)) {
                init_incid <- init_states[stem_object$dynamics$incidence_sources + 1]
                names(init_incid) <- names(stem_object$dynamics$incidence_codes)
                init_states <- c(init_states, init_incid)
                ode_incidence_codes <- stem_object$dynamics$incidence_codes + 1 # add 1 b/c of time column
        } else {
                ode_incidence_codes <- 0
        }

        ode_path <- stem_ode_path(ode_times        = ode_times,
                                  census_times     = census_times,
                                  census_inds      = census_inds,
                                  ode_pars         = ode_pars,
                                  init_states      = init_states,
                                  incidence_codes  = ode_incidence_codes,
                                  ode_pointer      = stem_object$dynamics$ode_pointers$ode_ptr,
                                  set_pars_pointer = stem_object$dynamics$ode_pointers$set_ode_params_ptr)

        colnames(ode_path) <- c("time",colnames(stem_object$dynamics$flow_matrix))

        return(ode_path)
}