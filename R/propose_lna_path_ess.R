#' Propose a latent LNA path with zero mean for use in the elliptical slice
#' sampler and return the paths of the LNA ODEs.
#'
#' @param path_cur list containing the current path and LNA ODEs
#' @param parameters list of lna parameters, the lna function pointer, and the
#'   stoichiometry matrix
#' @param times vector of times at which the LNA is to be evaluated
#' @param initdist_parameters vector of initial state parameters
#' @param fixed_inits is the initial state fixed or random?
#' @param fixed_inits are the initial states fixed?
#' @param initdist_parameters mean vector for the initial compartment counts
#' @param param_update_inds logical vector indicating at which of the times the
#'   LNA parameters need to be updated.
#' @param stoich_matrix the stoichiometry matrix, i.e. the transpose of the flow
#'   matrix.
#' @param drift_inds vector of indices corresponding to the drift process within
#'   the vector holding the ODEs
#' @param resid_inds vector of indices corresponding to the residual process
#'   within the vector holding the ODEs
#' @param diff_inds vector of indices corresponding to the diffusion process
#'   within the vector holding the ODEs
#' @param lna_pointer external pointer to the LNA function
#' @param log_scale logical
#' @param incidence_codes vector of incidence codes
#'
#' @return list with the proposed LNA path, along with the drift, residual, and
#'   diffusion processes
#' @export
propose_lna_path_ess <- function(path_cur, parameters, stoich_matrix, lna_ess_pointer, times, initdist_parameters, fixed_inits, param_update_inds, drift_inds, resid_inds, log_scale, incidence_codes = NULL) {

        # get number of times, compartments, and ODEs
        n_comps <- length(drift_inds)
        n_odes  <- 2 * n_comps
        n_times <- length(times)

        if(!is.null(incidence_codes)) {
                incid_codes_step <- incidence_codes + 1
                incid_codes_path <- incidence_codes + 2
                n_prev  <- n_comps - length(incidence_codes)
        }

        # create the objects to store the path and the ODEs the drift process is
        # the same as the current path but still needs to be reintegrated to
        # obtain the new residual path the diffusion path is unchanged
        path              <- matrix(0.0, nrow = n_times, ncol = n_comps + 1)    # path matrix
        path[,1]          <- as.numeric(times)
        drift_process     <- path_cur$path$drift_process                        # matrix with drift
        residual_process  <- matrix(0.0, nrow = n_times, ncol = n_comps)        # matrix for the residual process
        diffusion_process <- path_cur$path$diffusion_process                    # array with the diffusion matrices

        # initialize the path, drift process, and lna_state_vec
        if(fixed_inits) {
                path[1,-1]        <- drift_process[1,,drop = FALSE] # the path and drift are restarted at same point
        } else {
                stop("Non-fixed initial state not yet implemented.")
        }

        # set the lna_state vector
        lna_state_vec     <- c(drift_process[1,], residual_process[1,])

        # initialize the parameter vector
        lna_params <- list(parameters = parameters[1,], lna_ess_pointer = lna_ess_pointer, stoich = stoich_matrix)

        # iterate over the time sequence, solving the LNA over each interval
        for(j in 2:n_times) {

                # set the times of the interval endpoints
                t_L = times[j-1]
                t_R = times[j]

                # set the values of the
                if(param_update_inds[j]) lna_params$parameters = as.numeric(lna_pars[j-1,]);

                # integrate the LNA
                lna_state_vec <- deSolve::lsoda(y = lna_state_vec,
                                                func = CALL_LNA_ESS,
                                                times = c(t_L, t_R),
                                                parms = lna_params,
                                                maxsteps = 5e4)[2,-1]

                # transfer the elements of the lna_state_vec to the drift,
                # residual, and diffusion process objects
                residual_process[j,] <- lna_state_vec[resid_inds]

                # sample the next value ensure that there is not negative
                # incidence and that there are no negative values if the process
                # is not on the log scale
                lna_step <- mvtnorm::rmvnorm(1, drift_process[j,] + residual_process[j,], diffusion_process[,,j])

                path[j,-1]                <- lna_step                             # insert the path into the path matrix
                lna_state_vec[resid_inds] <- lna_step - drift_process[j,]         # adjust the residual process
        }

        return(list(path = path,
                    residual_path = path[,-1] - drift_process,
                    drift_process = drift_process,
                    residual_process = residual_process,
                    diffusion_process = diffusion_process))
}