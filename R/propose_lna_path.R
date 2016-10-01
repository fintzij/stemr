#' Propose a latent path using the linear noise approximation and return the
#' paths of the LNA ODEs.
#'
#' @param parameters list of lna parameters, the lna function pointer, and the
#'   stoichiometry matrix
#' @param times vector of times at which the LNA is to be evaluated
#' @param initdist_parameters vector of parameters for the initial state
#' @param fixed_inits is the initial state fixed or random?
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
#' @param add_noise add some noise to the diagonal of the diffusion matrix,
#'   defaults to FALSE
#'
#' @return list with the proposed LNA path, along with the drift, residual, and
#'   diffusion processes
#' @export
propose_lna_path <- function(parameters, stoich_matrix, lna_pointer, times, initdist_parameters, fixed_inits, param_update_inds, drift_inds, resid_inds, diff_inds, log_scale, incidence_codes = NULL) {

        # get number of times, compartments, and ODEs
        n_comps <- length(drift_inds)
        n_odes  <- diff_inds[length(diff_inds)]
        n_times <- length(times)

        if(!is.null(incidence_codes)) {
                incid_codes_step <- incidence_codes + 1
                incid_codes_path <- incidence_codes + 2
                n_prev  <- n_comps - length(incidence_codes)
        }

        # create the objects to store the path and the ODEs
        path              <- matrix(0.0, nrow = n_times, ncol = n_comps + 1)    # path matrix
        path[,1]          <- as.numeric(times)
        drift_process     <- matrix(0.0, nrow = n_times, ncol = n_comps)        # matrix for the drift process
        residual_process  <- matrix(0.0, nrow = n_times, ncol = n_comps)        # matrix for the residual process
        diffusion_process <- array(0.0, dim = c(n_comps, n_comps, n_times))     # array with the diffusion matrices

        # initialize the path, drift process, and lna_state_vec
        if(fixed_inits) {
                path[1,-1]        <- initdist_parameters
                drift_process[1,] <- path[1,-1]
                residual_path[1,] <- path[1,-1] - drift_process[1,]
        } else {
                stop("Random initial state not yet implemented.")
        }

        # set the lna_state vector
        lna_state_vec     <- c(drift_process[1,], numeric(length = n_comps + n_comps^2))

        # initialize the parameter vector
        lna_params <- list(parameters = parameters[1,], lna_pointer = lna_pointer, stoich = stoich_matrix)

        # iterate over the time sequence, solving the LNA over each interval
        for(j in 2:n_times) {

                # set the times of the interval endpoints
                t_L = times[j-1]
                t_R = times[j]

                # set the values of the
                if(param_update_inds[j]) lna_params$parameters = as.numeric(parameters[j-1,]);

                # integrate the LNA
                lna_state_vec <- deSolve::lsoda(y = lna_state_vec,
                                                func = CALL_COMPUTE_LNA,
                                                times = c(t_L, t_R),
                                                parms = lna_params,
                                                maxsteps = 5e4)[2,-1]

                # transfer the elements of the lna_state_vec to the drift,
                # residual, and diffusion process objects
                drift_process[j,]      <- lna_state_vec[drift_inds]
                residual_process[j,]   <- lna_state_vec[resid_inds]
                diffusion_process[,,j] <- lna_state_vec[diff_inds] + diag(1e-8, n_comps)

                # sample the next value ensure that there is not negative
                # incidence and that there are no negative values if the process
                # is not on the log scale
                if(log_scale) {
                        if(!is.null(incidence_codes)) {
                                lna_step <- try(tmvtnorm::rtmvnorm(1, drift_process[j,] + residual_process[j,], diffusion_process[,,j],
                                                                   lower = c(rep(-Inf, n_prev), path[j-1, incid_codes_path]),
                                                                   algorithm = "rejection"), silent = TRUE)
                                if(class(lna_step) == "try-error") {
                                        diffusion_process[,,j] <- as.matrix(Matrix::nearPD(diffusion_process[,,j],
                                                                                           ensureSymmetry = TRUE)$mat)
                                        lna_step <- tmvtnorm::rtmvnorm(1, drift_process[j,] + residual_process[j,],
                                                                       diffusion_process[,,j],
                                                                       lower = c(rep(-Inf, n_prev), path[j-1, incid_codes_path]),
                                                                       algorithm = "rejection")
                                }
                        } else {
                                lna_step  <- try(tmvtnorm::rtmvnorm(1, drift_process[j,] + residual_process[j,], diffusion_process[,,j],
                                                                algorithm = "rejection"), silent = TRUE)
                                if(class(lna_step) == "try-error") {
                                        diffusion_process[,,j] <- as.matrix(Matrix::nearPD(diffusion_process[,,j],
                                                                                          ensureSymmetry = TRUE)$mat)
                                        lna_step <- tmvtnorm::rtmvnorm(1, drift_process[j,] + residual_process[j,],
                                                                       diffusion_process[,,j],
                                                                       algorithm = "rejection")
                                        }
                        }
                } else if(!log_scale) {
                        if(!is.null(incidence_codes)) {
                                lna_step  <- try(tmvtnorm::rtmvnorm(1, drift_process[j,] + residual_process[j,], diffusion_process[,,j],
                                                                lower = c(rep(0, n_prev), path[j-1, incid_codes_path]),
                                                                algorithm = "rejection"), silent = TRUE)
                                if(class(lna_step) == "try-error") {
                                        diffusion_process[,,j] <- as.matrix(Matrix::nearPD(diffusion_process[,,j],
                                                                                           ensureSymmetry = TRUE)$mat)
                                        lna_step <- tmvtnorm::rtmvnorm(1, drift_process[j,] + residual_process[j,],
                                                                       diffusion_process[,,j],
                                                                       lower = c(rep(0, n_prev), path[j-1, incid_codes_path]),
                                                                       algorithm = "rejection")

                                }
                        } else {
                                lna_step  <- try(tmvtnorm::rtmvnorm(1, drift_process[j,] + residual_process[j,], diffusion_process[,,j],
                                                                lower = rep(0, n_comps),algorithm = "rejection"), silent = TRUE)
                                if(class(lna_step) == "try-error") {
                                        diffusion_process[,,j] <- as.matrix(Matrix::nearPD(diffusion_process[,,j],
                                                                                           ensureSymmetry = TRUE)$mat)
                                        lna_step <- tmvtnorm::rtmvnorm(1, drift_process[j,] + residual_process[j,],
                                                                       diffusion_process[,,j],
                                                                       lower = rep(0, n_comps),algorithm = "rejection")
                                }
                        }
                }

                path[j,-1]                <- lna_step                             # insert the path into the path matrix
                lna_state_vec[resid_inds] <- lna_step - lna_state_vec[drift_inds] # set the residual term
                lna_state_vec[diff_inds]  <- 0.0                                  # zero out the diffusion process
        }

        return(list(path = path,
                    residual_path = path[,-1] - drift_process,
                    drift_process = drift_process,
                    residual_process = residual_process,
                    diffusion_process = diffusion_process))
}