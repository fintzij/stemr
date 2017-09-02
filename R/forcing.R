#' Declare a time varying covariate to be a forcing variable that moves
#' individuals between model compartments at discrete times.
#'
#' @param tcovar_name name of the time varying covariate
#' @param from name of the compartment that individuals exit
#' @param to name of the compartment that individuals enter
#'
#' @return list containing the forcing specification
#' @export
#'
#' @examples forcings = list(forcing("Su2Sv", "S_unvaccinated", "S_vaccinated"))
forcing <- function(tcovar_name, from, to) {
        list(tcovar_name = tcovar_name, from = from, to = to)
}