#' Declare a time varying covariate to be a forcing variable that moves
#' individuals between model compartments at discrete times. Flow in and out of
#' model compartments is allocated proportionally to the compartment counts in
#' the source and destination compartments.
#'
#' @param tcovar_name name of the time varying covariate
#' @param from vector of compartment names from which individuals exit
#' @param to vector of compartment names to which individuals enter
#'
#' @return list containing the forcing specification
#' @export
forcing <- function(tcovar_name, from, to) {
      
      if(length(from) != length(to)) stop("The number of forcings out must equal the number of forcings in.")
      
      list(
            tcovar_name = tcovar_name,
            from = from,
            to = to
      )
}