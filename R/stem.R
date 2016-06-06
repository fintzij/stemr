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
#' @param parameters Named vector of parameter values.
#' @param strata names of model strata.
#' @param constants Named vector of constants.
#' @param stem_settings otional list of inference settings, most
#'   straightforwardly generated using the \code{stem_control} function.
#' @param meas_process list of functions to simulate from or evaluate the
#'   likelihood of the measurement process. These are most easily generated
#'   using the \code{stem_measure} function.
#'
#' @return returns a \code{stem} list.
#' @export
#'
stem <- function(stem_object = NULL, data = NULL, dynamics = NULL, stem_settings = NULL, meas_process = NULL) {

        # if a dataset or multiple datasets were provided, collect the
        # observation times and indicator matrices for which compartments were
        # measured at which times
        if(!is.null(data)) {
                if(!is.list(data)) {
                        times <- data[, 1, drop = FALSE]
                        comps_at_obstimes <- matrix(1, nrow = nrow(data), ncol = ncol(data) - 1)
                        colnames(comps_at_obstimes) <- colnames(data)[-1]
                }
                else{
                        dat_times <- lapply(data, function(x) x[,1]); times <- as.matrix(sort(unlist(dat_times)))
                        comps_at_obstimes <- matrix(0, nrow = length(times), ncol = sum(sapply(data, function(x) ncol(x)-1)))
                        colnames(comps_at_obstimes) <- sapply(data, function(x) colnames(x)[-1])
                        for(s in 1:ncol(comps_at_obstimes)) {
                                comps_at_obstimes[times %in% dat_times[[s]], s] <- 1
                        }
                }
        }

        if(is.null(stem_object)) {
                stem_object <- structure(list(data              = NULL,
                                              times             = NULL,
                                              comps_at_obstimes = NULL,
                                              dynamics          = NULL,
                                              meas_process      = NULL,
                                              stem_settings     = NULL), class = "stem")
        }

        if(!is.null(data))              stem_object$data              <- data
        if(!is.null(times))             stem_object$times             <- times
        if(!is.null(comps_at_obstimes)) stem_object$comps_at_obstimes <- comps_at_obstimes
        if(!is.null(dynamics))          stem_object$dynamics          <- dynamics
        if(!is.null(stem_settings))     stem_object$stem_settings     <- stem_settings

        return(stem_object)
}