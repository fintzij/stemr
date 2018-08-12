#' Parses a function of each of a vector of model compartments and aggregates
#' the results.
#'
#' This function serves as a placeholder in order to generate documentation on
#' how to specify the comp_fcn string. It also checks that the aggregation term
#' is correctly specified.
#'
#' \code{comp_fcn} facilitates quick specification of rates that depend on
#' multiple compartments. It is provided as part of a rate string, and is parsed
#' and evaluated internally in order to generate a usable function of model
#' compartments. In particular, it allows a user to consisely specify a function
#' that is to be applied to each desired model compartment, and whether the
#' resulting quantities should be added or multiplied together.
#'
#' @param fcn function of each the model compartments, e.g. log(I_ADJ) would
#'   apply the log function to each of the adjacent compartments prior to
#'   aggregation.
#' @param aggregation either sum or prod for whether the sum or product of the
#'   model compartments is desired.
#' @param compartments list of compartments to be substituted for the reserved
#'   word, e.g. all adjacent compartments would be substituted for I_ADJ. This
#'   argument is generated internally and should not be provided by the user.
#'
#' @return Nothing. Text string is parsed internally by \code{sub_comp_rate} and
#'   \code{sub_comp_emit}.
#' @export
#'
#' @examples # The following rate specification would be used to model the force
#' # of infection as the sum of beta times the number of infecteds in an
#' # individual's own compartment, and gamma times the sum of the log number of
#' # infecteds in adjacent compartments.
comp_fcn <- function(fcn, aggregation, compartments = NULL) {

        if(!aggregation %in% c("sum", "prod")) {
                stop("The aggregation for model compartments must be given as either sum or prod.")
        }

}