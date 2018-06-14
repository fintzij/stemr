#' Define a parameter block for an MCMC kernel
#'
#' @param pars either a numeric vector of parameter indices or a character
#'   vector of parameter names on their estimation scale (i.e., corresponding to
#'   row/column names in the kernel covariance matrix)
#'
#' @return parameter block for use in MCMC kernel
#' @export
parblock <- function(pars) {
      return(list(pars = pars))
}