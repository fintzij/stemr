#' Simulate a path using the linear noise approximation
#'
#' @param stem_object stochastic epidemic model object
#' @inheritParams simulate_stem
#'
#' @return matrix containing a path simulated using the LNA
#' @export
simulate_lna <- function(stem_object, census_times, lna_restart, init_states) {

        # initialize the path matrix and the deterministic, stochastic, and
        # residual processes
        lna_path        <- matrix(0, nrow = length(census_times), ncol = 1 + length(init_states)) # create matrix
        lna_path[1,]    <- c(census_times[1], init_states)        # assign initial state values
        lna_path[,1]    <- census_times                           # insert the census times

        # pull the incidence codes and compartment codes
        comp_codes <- stem_object$dynamics$comp_codes + 1

        if(is.null(stem_object$dynamics$incidence_codes)) {
                incidence_codes <- NULL

        } else {
                incidence_codes   <- stem_object$dynamics$incidence_codes + 1
        }

        # if the restarting version is to be used, grab the restart times
        if(lna_restart) {
                if(is.logical(lna_restart)) {
                        restart_times <- census_times
                } else {
                        restart_times <- c(lna_restart, census_times[length(census_times)])
                        if(restart_times[1] != census_times[1]) restart_times <- c(census_times[1], restart_times)
                }

                r <- 2
                r_time <- restart_times[r]
        }

        # initialize the deterministic, stochastic, and residual processes
        det_proc        <- init_states                                          # deterministic process
        innov_proc      <- det_proc * 0                                         # stochastic process
        resid_proc      <- matrix(0, length(det_proc), length(det_proc))        # residual process
        n_comps         <- length(stem_object$dynamics$comp_codes) + length(stem_object$dynamics$incidence_codes) # number of model compartments

        # initialize the model list
        stemr_lnamod    <- list(parameters  = stem_object$dynamics$parameters,
                                constants   = stem_object$dynamics$constants,
                                tcovar      = stem_object$dynamics$lna_tcovar[1, -1, drop = FALSE],
                                flow_matrix = t(stem_object$dynamics$flow_matrix),
                                hazard_ptr  = stem_object$dynamics$lna_ptrs$hazard_ptr,
                                jacobian_ptr = stem_object$dynamics$lna_ptrs$jacobian_ptr)

        for(k in 2:(length(census_times))) {

                # increment the interval endpoints
                t_L <- census_times[k-1]
                t_R <- census_times[k]

                # retrieve the time-varying covariate values
                stemr_lnamod$tcovar <- stem_object$dynamics$lna_tcovar[k-1, -1, drop = FALSE]

                # solve the lna
                lna_soln <- solve_lna(stemr_lnamod, det_proc, innov_proc, t_L, t_R, n_comps)

                # if there are incidence compartments, the noise needs to be simulated just based on the non-incidence variables
                # simulate the next value - redraw if any components are negative
                lna_step <- MASS::mvrnorm(1, lna_soln$det_proc + lna_soln$innov_proc, lna_soln$resid_proc)

                # while(any(lna_step < 0) || (!is.null(incidence_codes) && any(lna_step[incidence_codes] < lna_path[k-1,incidence_codes + 1]))) {
                while(any(lna_step < 0)) {
                        lna_step <- MASS::mvrnorm(1, lna_soln$det_proc + lna_soln$innov_proc, lna_soln$resid_proc)
                }

                # insert the new path into the path matrix
                lna_path[k,-1] <- lna_step

                # update the deterministic process, the innovation, and the residual process
                if(!lna_restart) {
                        det_proc   <- lna_soln$det_proc
                        innov_proc <- lna_step - det_proc

                } else if(t_R == r_time){
                        # restart the LNA
                        det_proc   <- lna_step
                        innov_proc <- rep(0, n_comps)

                        # increment to the next restart time
                        r          <- r+1
                        r_time     <- restart_times[r]

                } else if(t_R != r_time) {
                        det_proc   <- lna_soln$det_proc
                        innov_proc <- lna_step - det_proc
                }
        }

        colnames(lna_path) <- c("time", names(stem_object$dynamics$comp_codes), names(stem_object$dynamics$incidence_codes))

        # if some of the compartments track incidence, compute the incidence
        if(!is.null(stem_object$dynamics$incidence_codes)) {
                census_incidence_rows <- rep(list(seq_along(census_times) - 1), length(stem_object$dynamics$incidence_codes))
                compute_incidence(censusmat = lna_path, col_inds  = incidence_codes, row_inds  = census_incidence_rows)
                for(i in incidence_codes) lna_path[,i + 1] <- pmax.int(lna_path[,i + 1], 0)
        }

        return(lna_path)
}