#' Solve the LNA ODEs over a single time interval
#'
#' @param stemr_lnamod list stemr LNA model objects
#' @param det_proc vector for the deterministic process
#' @param innov_proc vector for the stochastic term
#' @param t_L left endpoint of the interval
#' @param t_R right endpoint of the interval
#' @param n_comps number of model compartments
#'
#' @return list containing the solution to the LNA over the time interval
#'   [t_L,t_R].
#' @export
solve_lna <- function(stemr_lnamod, det_proc, innov_proc, t_L, t_R, n_comps) {

        lna_soln <- deSolve::lsoda(y = c(det_proc, innov_proc, rep(0, n_comps^2)),
                                   func = lna_odes,
                                   times = c(t_L, t_R),
                                   parms = stemr_lnamod)[2,-1]

        return(list(det_proc = lna_soln[1:n_comps],
                    innov_proc = lna_soln[(n_comps + 1):(2*n_comps)],
                    resid_proc = matrix(lna_soln[(2*n_comps + 1):length(lna_soln)], ncol = n_comps)))
}