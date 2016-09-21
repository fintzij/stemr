#' Compute the log density of an LNA path by integrating the LNA ODEs.
#'
#' @param path matrix containing a sampled LNA path
#' @inheritParams propose_lna_path
#' @param fixed_inits is the initial state fixed
#'
#' @return log-density of the LNA path
#' @export
compute_lna_density2 <- function(path, parameters, stoich_matrix, lna_pointer, times, param_update_inds, drift_inds, resid_inds, diff_inds, log_scale, incidence_codes, fixed_inits) {

        # get number of times, compartments, and ODEs
        n_comps <- length(drift_inds)
        n_odes  <- diff_inds[length(diff_inds)]
        n_times <- length(times)

        if(!is.null(incidence_codes)) {
                incid_codes_step <- incidence_codes + 1
                incid_codes_path <- incidence_codes + 2
                n_prev  <- n_comps - length(incidence_codes)
        }

        # create the objects to store the ODEs
        drift_process     <- matrix(0.0, nrow = n_times, ncol = n_comps)        # matrix for the drift process
        residual_process  <- matrix(0.0, nrow = n_times, ncol = n_comps)        # matrix for the residual process
        diffusion_process <- array(0.0, dim = c(n_comps, n_comps, n_times))     # array with the diffusion matrices

        # initialize the path, drift process, and lna_state_vec
        drift_process[1,] <- path$path[1,-1]

        # initialize the log likelihood
        if(fixed_inits) {
                log_lik <- 0
        } else {
                log_lik <- tmvtnorm::dtmvnorm(x = path$path[1,-1],
                                              mean = path$drift_process[1,],
                                              sigma = path$diffusion_process[,,1],
                                              lower = ifelse(log_scale, rep(-Inf, n_comps), rep(0, n_comps)),
                                              log = T)
        }

        # set the lna_state vector
        lna_state_vec  <- c(drift_process[1,], numeric(length = n_comps + n_comps^2))

        # initialize the parameter vector
        lna_params     <- list(parameters = parameters[1,], lna_pointer = lna_pointer, stoich = stoich_matrix)

        # iterate over the time sequence, solving the LNA over each interval
        for(j in 2:n_times) {

                # set the times of the interval endpoints
                t_L = times[j-1]
                t_R = times[j]

                # set the values of the
                if(param_update_inds[j]) lna_params$parameters = as.numeric(lna_pars[j-1,]);

                # integrate the LNA
                lna_state_vec <- deSolve::lsoda(y = lna_state_vec,
                                                func = CALL_COMPUTE_LNA,
                                                times = c(t_L, t_R),
                                                parms = lna_params,
                                                maxsteps = 5e4)[2,-1]

                # transfer the elements of the lna_state_vec to the drift,
                # residual, and diffusion process objects
                vec2procmats(lna_state_vec, drift_process, residual_process, diffusion_process, j)

                if(!all(drift_process[j,] == path$drift_process[j,])) break

                # sample the next value ensure that there is not negative
                # incidence and that there are no negative values if the process
                # is not on the log scale
                if(log_scale) {
                        if(!is.null(incidence_codes)) {
                                log_lik <- log_lik + tmvtnorm::dtmvnorm(path$path[j, -1],
                                                                      drift_process[j,] + residual_process[j,],
                                                                      diffusion_process[,,j],
                                                                      lower = c(rep(-Inf, n_prev), path$path[j-1, incid_codes_path]),
                                                                      log = TRUE)
                        } else {
                                log_lik <- log_lik + tmvtnorm::dtmvnorm(path$path[j, -1],
                                                                      drift_process[j,] + residual_process[j,],
                                                                      diffusion_process[,,j],
                                                                      log = TRUE)
                        }

                } else if(!log_scale) {
                        if(!is.null(incidence_codes)) {
                                log_lik <- log_lik + tmvtnorm::dtmvnorm(path$path[j, -1],
                                                                      drift_process[j,] + residual_process[j,],
                                                                      diffusion_process[,,j],
                                                                      lower = c(rep(0, n_prev), path$path[j-1, incid_codes_path]),
                                                                      log = TRUE)
                        } else {
                                log_lik <- log_lik + tmvtnorm::dtmvnorm(path$path[j, -1],
                                                                      drift_process[j,] + residual_process[j,],
                                                                      diffusion_process[,,j],
                                                                      lower = rep(0, n_comps),
                                                                      log = TRUE)
                        }
                }

                lna_state_vec[resid_inds] <- path$path[j, -1] - lna_state_vec[drift_inds] # set the residual term
                lna_state_vec[diff_inds]  <- 0.0                                  # zero out the diffusion process
        }

        return(log_lik)
}