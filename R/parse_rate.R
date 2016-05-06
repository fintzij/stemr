#' Transform a text snippet for a rate function into C++ text for the rate
#' function.
#'
#' @param rate string snippet declaring how to compute the rate, must be a
#'   function of parameters, compartment_strata, and constants, and must be
#'   valid C++ code.
#' @param parameters vector of parameters
#' @param compartment_names vector of compartment names
#' @param constants vector of constant names
#' @param strata character vector of strata names to which the rate applies, or
#'   numeric vector of strata numbers
#' @param adjacency adjacency matrix
#'
#' @return string that can be evaluated and compiled to create a C++ function
#'   for computing the rate.
#' @export
#'
parse_rate <- function(rate, parameters, compartment_names, constants, strata, adjacency) {



        return(rate_function)
}