#' Assign a vector of parameters to a stem object
#'
#' @param parameters named vector of parameters 
#' @param stem_object stem list
#'
#' @return stem_object with new parameters
#' @export
stem_parameters <- 
      function(parameters, stem_object) {
            stem_object$dynamics$parameters = parameters
            return(stem_object)
}