#' Construct a stem object.
#'
#' @param data,times The observation times and data. The user may specify either
#'   a vector of observation times in the \code{times} argument, or a list of
#'   vectors of observation times, where the list element names identify the
#'   compartment_strata to which the observation times correspond (e.g.
#'   I_children, I_adults). \code{times} must be numeric and strictly
#'   increasing.
#'
#'   \code{data} may be provided either as a data frame of observations with
#'   number of rows  equal to the number of observation times and with column
#'   names that identify the compartment_stratum (e.g. I_children, I_adults)
#'   being measured, or as a list of vectors, where the list element names
#'   identify the compartment_stratum and where each vector has a corresponding
#'   vector of observation times.
#' @param dynamics A list of objects describing the model dynamics, most
#'   straighforwardly generated using the \code{stem_dynamics} function.
#' @param parameters Named vector of parameter values.
#' @param constants Names vector of constants. If constants are stratum
#'   specific, they must be named appropriately - e.g. N_children, N_adults.
#' @param strata Names of strata. Must be specified if there are multiple
#'   strata.
#' @param timevar name of the time variable
#' @param stem_settings otional list of simulation settings, most
#'   straightforwardly generated using the \code{stem_control} function.
#' @param rmeasure,dmeasure functions to simulate from or evaluate the
#'   likelihood of the measurement process. At least one must be specified.
#' @param incidence_vars vector of compartments for which the data consists of
#'   incidence data.
#'
#' @return returns a \code{stem} list.
#' @export
#'
stem <- function(data = NULL, times = NULL, dynamics, parameters, constants = NULL, strata = NULL, timevar = NULL, stem_settings = NULL, rmeasure = NULL, dmeasure = NULL, incidence_vars = NULL) {


}