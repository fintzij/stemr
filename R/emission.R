#' Generates an emission list to be supplied to the \code{\link{stem_measure}}
#' function.
#'
#' @param compartment compartment being measured or simulated from, a string.
#' @param obstimes vector of observation times.
#' @param distribution emission distribution, one of "binomial", "poisson",
#'   "negbinomial", "gaussian"
#' @param parameters parameter names of the emission distribution, a string
#'   vector
#' @param strata
#'
#' @return
#' @export
#'
#' @examples
emission <- function(compartment, obstimes, distribution, parameters, strata = NULL) {

        # generate the list
        list(compartment = compartment, obstimes = obstimes, distribution = distribution, parameters = parameters, strata = strata)
}