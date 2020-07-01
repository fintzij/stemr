#' Set the parameter values for a stochastic epidemic model object.
#'
#' @param stem_object stochastic epidemic model object for which parameter
#'   values should be set.
#' @param ...
#'
#' @return stem object with updated parameters
#' @export
set_params <- function(stem_object, ...) {

        newpars      <- list(...)
        newpar_names <- lapply(newpars, names)
        oldpar_names <- names(stem_object$dynamics$parameters)

        for(l in seq_along(newpars)) {
                stem_object$dynamics$parameters[match(newpar_names[[l]], oldpar_names)] <- newpars[[l]]
        }

        return(stem_object)
}