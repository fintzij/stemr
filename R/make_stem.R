#' Construct a stem object.
#'
#' @param stem_object an existing stochastic epidemic model object to which
#'   model settings should be added.
#' @param data a matrix/data frame, or a list of matrices/data frames. All
#'   columns must be named according to which compartment_strata are measured.
#'   The first column must consist of observation times, t_1,...,t_L. If data on
#'   all measured compartments are accrued at the same observation times, a
#'   single matrix or data frame may be provided. If each compartment in the
#'   data was measured at different observation times, a list of matrices or
#'   data frames must be provided. Again, the first column of each matrix must
#'   consist of observation times, while subsequent columns must be labeled
#'   according to which compartments are reflected.
#' @param dynamics A list of objects describing the model dynamics, most
#'   straighforwardly generated using the \code{stem_dynamics} function.
#' @param measurement_process list of functions to simulate from or evaluate the
#'   likelihood of the measurement process. These are most easily generated
#'   using the \code{stem_measure} function.
#'
#' @return returns a \code{stem} object.
#' @export
make_stem <- function(stem_object = NULL, data = NULL, dynamics = NULL, measurement_process = NULL) {

        gc(full = T)
        
        if(is.null(stem_object)) {
                stem_object <- 
                        list(dynamics            = NULL,
                             measurement_process = NULL)
        }

        if(!is.null(dynamics)) {
                if(!is.null(stem_object$dynamics)) {
                        stem_object$dynamics = NULL
                }
                gc(full = T)
                stem_object$dynamics <- dynamics
        }
        
        if(!is.null(measurement_process)) {
                if(!is.null(stem_object$measurement_process)) {
                        stem_object$measurement_process = NULL
                }
                gc(full = T)
                stem_object$measurement_process <- measurement_process
        } 
        
        return(stem_object)
}
