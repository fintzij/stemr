#' Simulate a path using the linear noise approximation
#'
#' @param nsim integer for the number of simulations
#' @param lna_times vector of times at which the LNA must be evaluated
#' @param census_times vector of times at which the LNA should be censused
#' @param census_inds logical vector the same length as lna_times indicating
#'   whether the LNA path is censused at each of the lna times
#' @param lna_pars numeric matrix of parameters, constants, and time-varying
#'   covariates at each of the lna_times
#' @param init_states matrix of initial state vectors
#' @param incidence_codes vector of incidence codes, or 0 if no incidence
#'   compartments
#' @param log_scale is the LNA applied on the log scale?
#' @param flow_matrix flow matrix for the stochastic epidemic model object
#' @param lna_pointer external pointer to LNA integration fcn
#'
#' @return matrix (array if nsim>1) containing a path(s) simulated using the LNA
#' @export
simulate_lna <- function(nsim, lna_times, census_times, census_inds, lna_pars, init_states, incidence_codes, log_scale, flow_matrix, lna_pointer) {

        # get the dimensions of various objects
        n_comps <- ncol(init_states)           # number of model compartments
        n_odes  <- 2*n_comps + n_comps*n_comps # number of ODEs
        n_times <- length(lna_times)           # number of times at which the LNA must be evaluated
        n_incidence <- length(incidence_codes) # number of incidence codes
        n_census_times <- length(census_times) # number of census times

        # initialize the objects used in each time interval
        t_L <- lna_times[1]
        t_R <- lna_times[2]
        lna_params <- list(parameters = lna_pars[1,], lna_pointer = lna_pointer, stoich = t(flow_matrix))

        # create an armadillo cube to store the simulations
        lna_paths <- array(0, dim = c(n_census_times, n_comps + 1, nsim))
        lna_paths[,1,] <- census_times

        # full_lna_paths <- array(0, dim = c(n_census_times, n_comps*2 + n_comps^2 + 1, nsim))

        # initialize the LNA objects - the vector for storing the ODES, the state vector, and the Jacobian
        lna_state_vec     <- numeric(length = n_odes)                # vector to store the results of the ODEs
        drift_process     <- numeric(length = n_comps)               # vector for the deterministic drift process
        residual_process  <- numeric(length = n_comps)               # vector for the residual process
        diffusion_process <- matrix(as.numeric(0), n_comps, n_comps) # matrix for the diffusion

        # indices in the ODE vector for the drift, residual, and diffusion
        drift_inds <- seq_len(n_comps);
        resid_inds <- seq_len(n_comps) + n_comps;
        diff_inds  <- seq_len(n_comps * n_comps) + 2*n_comps;

        # start simulating the LNA paths
        for(k in 1:nsim) {

                # set the initial values of the LNA
                lna_paths[1,-1,k]    <- as.numeric(init_states[k,]); # assign the initial state to the path matrix
                drift_process        <- as.numeric(init_states[k,]); # assign the initial state to the drift vector
                residual_process[]   <- 0.0;               # zero out the residual vector
                diffusion_process[,] <- 0.0;               # zero out the diffusion matrix

                procs2vec(lna_state_vec, drift_process, residual_process, diffusion_process)

                # full_lna_paths[1,,k] <- c(t_L,lna_state_vec)

                # iterate over the time sequence, solving the LNA over each interval
                for(j in 2:n_times) {

                        # set the times of the interval endpoints
                        t_L = lna_times[j-1];
                        t_R = lna_times[j];

                        # set the values of the
                        lna_params$parameters = as.numeric(lna_pars[j-1,]);

                        # integrate the LNA
                        lna_state_vec <- deSolve::lsoda(y = lna_state_vec,
                                                       func = CALL_COMPUTE_LNA,
                                                       times = c(t_L, t_R),
                                                       parms = lna_params,
                                                       maxsteps = 5e4)[2,-1]

                        # transfer the elements of the lna_state_vec to the process objects
                        vec2procs(lna_state_vec, drift_process, residual_process, diffusion_process)

                        # sample the next value
                        # lna_step <- drift_process + residual_process
                        lna_step <- mvtnorm::rmvnorm(1, drift_process + residual_process, diffusion_process, method = "svd")
                        # lna_step <- mvtnorm::rmvnorm(1, drift_process + residual_process, diffusion_process,method = "svd")

                        # if any values are negative and the process is on the natural scale, resample them
                        while(!log_scale && any(lna_step < 0)) {
                                lna_step <- mvtnorm::rmvnorm(1, drift_process + residual_process, diffusion_process, method = "svd")
                                # lna_step <- mvtnorm::rmvnorm(1, drift_process + residual_process, diffusion_process,method = "svd")
                        }

                        # insert the sampled value into the path if called for
                        if(census_inds[j]) {
                                lna_paths[j,-1,k] <- lna_step
                        }

                        # restart the LNA if called for
                        # if(restart_inds[j]) {
                        #         drift_process        <- lna_step[1,]   # restart the deterministic process at the sampled state
                        #         # drift_process        <- lna_step   # restart the deterministic process at the sampled state
                        #         # full_lna_paths[j,,k] <- c(t_R,lna_state_vec)
                        #         residual_process[]   <- 0.0            # set the residual process to zero
                        #         diffusion_process[,] <- 0.0            # the diffusion process is set to zero
                        # } else {
                                # the residual process is the difference between the deterministic
                                # process and the sampled state and the diffusion process is set to zero
                                residual_process     <- lna_step - drift_process
                                # residual_process     <- lna_step - drift_process
                                # full_lna_paths[j,,k] <- c(t_R,lna_state_vec)

                                diffusion_process[,] <- 0.0
                        # }

                        # copy the process objects back into the lna state vector
                        procs2vec(lna_state_vec, drift_process, residual_process, diffusion_process);
                }
        }

        # compute the incidence if necessary
        if(incidence_codes[1] != 0) {

                # iterate through the path array
                for(k in 1:nsim) {
                        lna_paths[1,incidence_codes,k] <- 0
                        lna_paths[2:n_census_times, incidence_codes, k] <- diff(lna_paths[, incidence_codes, k])

                        # clamp the slice to make sure there are no negative values if the process is on the linear scale
                        if(!log_scale) {
                                lna_paths[,,k] = pmax(lna_paths[,,k], 0.0);
                        }
                }
        }

        return(lna_paths)
}