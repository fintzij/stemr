#' Generate a list to be used in specifying a time-varying parameter that has a 
#' latent Gaussian distribution and is updated via elliptical slice sampling.
#' 
#' @param tparam_name name of the time--varying parameter
#' @param times vector of times when the time-varying parameter changes.
#' @param draws2par function for mapping a vector of N(0,1) draws of length 
#'   equal to the length of the \code{times} argument. The function should take 
#'   two arguments, the vector of model hyperparameters (i.e. the parameters 
#'   argument) and a vector of N(0,1) draws, and return a vector of time-varying
#'   parameter values. Note that the function should define a deterministic 
#'   transformation.
#' @param values vector of values of N(0,1) draws for the time-varying 
#'   parameter, defaults to a vector of zeros. The values are computed by
#'   applying the \code{draws2par} function to the supplied vector. 
#'   
#' @return list to be used in specifying a time-varying parameter. \describe{A 
#'   time-varying parameter is defined as a (possibly non-linear) function of a 
#'   set of N(0,1) draws that are updated via ellipeical slice sampling.}
#' @export
tpar <- function(tparam_name, times, draws2par, values = NULL) {
      
      if(is.null(values)) values <- rep(0.0, length(times))
      
      return(list(tparam_name = tparam_name, times = times, draws2par = draws2par, values = values))
}