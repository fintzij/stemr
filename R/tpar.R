#' Generate a list to be used in specifying a time-varying parameter that has a
#' latent Gaussian distribution and is updated via elliptical slice sampling.
#'
#' @param tparam_name name of the time--varying parameter
#' @param draws2par function for mapping a vector of N(0,1) draws of length
#'   equal to the length of the \code{times} argument. The function should have
#'   the following signature: \code{draws2par(draws, pars)}. Here, draws is a
#'   vector of N(0,1) draws and pars is a vector of parameters. The
#'   function should return a vector of time-varying parameter values.
#' @param times vector of times when the time-varying parameter changes.
#' @param n_draws number of N(0,1) random variates
#' @param values vector of values of N(0,1) draws for the time-varying
#'   parameter, defaults to a vector of zeros. The values are computed by
#'   applying the \code{draws2par} function to the supplied vector.
#'
#' @return list to be used in specifying a time-varying parameter. \describe{A
#'   time-varying parameter is defined as a (possibly non-linear) function of a
#'   set of N(0,1) draws that are updated via elliptical slice sampling.}
#' @export
tpar <- function(tparam_name, draws2par, times, n_draws, values = NULL) {
      
      if(is.null(values)) values <- rep(0.0, length(times))
      
      if(is.null(control)) control = tpar_control()
      
      return(list(
          tparam_name = tparam_name,
          draws2par = draws2par,
          times = times,
          n_draws = n_draws,
          values = values
      ))
}