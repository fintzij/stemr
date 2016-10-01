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
        # get number of times, compartments, and ODEs
        n_comps    <- ncol(flow_matrix)
        n_odes     <- 2*n_comps + n_comps^2
        n_times    <- length(lna_times)
        n_census_times <- length(census_times)

        if(incidence_codes[1] != 0) {
                incid_codes_step <- incidence_codes + 1
                incid_codes_path <- incidence_codes + 2
                n_prev  <- n_comps - length(incidence_codes)
        }

        # create an armadillo cube to store the simulations
        lna_paths <- array(0, dim = c(n_census_times, n_comps + 1, nsim))
        lna_paths[,1,] <- census_times

        # initialize the objects used in each time interval
        t_L <- lna_times[1]
        t_R <- lna_times[2]
        lna_params <- list(parameters = lna_pars[1,], lna_pointer = lna_pointer, stoich = t(flow_matrix))

        # initialize the LNA objects - the vector for storing the ODES, the state vector, and the Jacobian
        lna_state_vec     <- numeric(length = n_odes)                # vector to store the results of the ODEs
        drift_process     <- numeric(length = n_comps)               # vector for the deterministic drift process
        residual_process  <- numeric(length = n_comps)               # vector for the residual process
        diffusion_process <- matrix(as.numeric(0), n_comps, n_comps) # matrix for the diffusion

        # start simulating the LNA paths
        for(k in 1:nsim) {

                # set the initial values of the LNA
                lna_paths[1,-1,k]    <- as.numeric(init_states[k,]); # assign the initial state to the path matrix
                drift_process        <- as.numeric(init_states[k,]); # assign the initial state to the drift vector
                residual_process[]   <- 0.0;               # zero out the residual vector
                diffusion_process[,] <- 0.0;               # zero out the diffusion matrix

                procs2vec(lna_state_vec, drift_process, residual_process, diffusion_process)

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

                        # enforce diagonal dominance
                        diffusion_process <- diffusion_process + diag(1e-8, n_comps)

                        # sample the next value
                        if(log_scale) {
                                if(incidence_codes[1] != 0) {
                                        lna_step <- try(tmvtnorm::rtmvnorm(1, drift_process + residual_process, diffusion_process,
                                                                           lower = c(rep(-Inf, n_prev),lna_paths[j-1,incid_codes_path,k]),
                                                                           algorithm = "rejection"), silent = TRUE)
                                        if(class(lna_step) == "try-error") {
                                                diffusion_process <-as.matrix(Matrix::nearPD(diffusion_process,ensureSymmetry = TRUE)$mat)
                                                lna_step <- tmvtnorm::rtmvnorm(1, drift_process + residual_process, diffusion_process,
                                                                               lower = c(rep(-Inf, n_prev),
                                                                                         lna_paths[j-1,incidence_codes_path,k]),
                                                                               algorithm = "rejection")
                                        }
                                } else {
                                        lna_step  <- try(tmvtnorm::rtmvnorm(1, drift_process + residual_process, diffusion_process,
                                                                            algorithm = "rejection"), silent = TRUE)
                                        if(class(lna_step) == "try-error") {
                                                diffusion_process <-as.matrix(Matrix::nearPD(diffusion_process,ensureSymmetry = TRUE)$mat)
                                                lna_step <- tmvtnorm::rtmvnorm(1, drift_process + residual_process, diffusion_process,
                                                                               algorithm = "rejection")
                                        }
                                }
                        } else if(!log_scale) {
                                if(incidence_codes[1] != 0) {
                                        lna_step  <- try(tmvtnorm::rtmvnorm(1, drift_process + residual_process, diffusion_process,
                                                                            lower = c(rep(0, n_prev), path[j-1, incid_codes_path, k]),
                                                                            algorithm = "rejection"), silent = TRUE)
                                        if(class(lna_step) == "try-error") {
                                                diffusion_process <-as.matrix(Matrix::nearPD(diffusion_process,ensureSymmetry = TRUE)$mat)
                                                lna_step <- tmvtnorm::rtmvnorm(1, drift_process + residual_process, diffusion_process,
                                                                               lower = c(rep(0, n_prev), path[j-1, incid_codes_path, k]),
                                                                               algorithm = "rejection")

                                        }
                                } else {
                                        lna_step  <- try(tmvtnorm::rtmvnorm(1, drift_process + residual_process, diffusion_process,
                                                                            lower = rep(0, n_comps),algorithm = "rejection"), silent = TRUE)
                                        if(class(lna_step) == "try-error") {
                                                diffusion_process <-as.matrix(Matrix::nearPD(diffusion_process,ensureSymmetry = TRUE)$mat)
                                                lna_step <- tmvtnorm::rtmvnorm(1, drift_process + residual_process, diffusion_process,
                                                                               lower = rep(0, n_comps),algorithm = "rejection")
                                        }
                                }
                        }

                        # insert the sampled value into the path if called for
                        if(census_inds[j]) {
                                lna_paths[j,-1,k] <- lna_step
                        }

                        residual_process     <- c(lna_step - drift_process)
                        diffusion_process[,] <- 0.0

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
                }
        }

        return(lna_paths)
}