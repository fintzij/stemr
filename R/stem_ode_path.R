#' Solve the ODEs over a single time interval
#'
#' @param stemr_odemod list stemr ODE model objects
#' @param ode_state vector for the deterministic process
#' @param t_L left endpoint of the interval
#' @param t_R right endpoint of the interval
#'
#' @return list containing the solution to the ODE over the time interval
#'   [t_L,t_R].
#' @export
stem_ode_path <- function(stemr_odemod, ode_state, t_L, t_R) {

        return(deSolve::lsoda(y = ode_state,
                              func = stemr_odemod$ode_fcn,
                              times = c(t_L, t_R),
                              parms = stemr_odemod)[2,-1])
}