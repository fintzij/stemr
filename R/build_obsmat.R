#' Generate a template for a stemr observation matrix, or combine multiple
#' datasets into a stemr observation matrix; used internally.
#'
#' @param meas_procs parsed list of emission lists, generated internally within
#'   \code{\link{stem_measure}}.
#' @param datasets list of datasets.
#'
#' @return observation matrix
#' @export
build_obsmat <- function(meas_procs = NULL, datasets = NULL) {

        # if a list of emission lists has been provided
        if(!is.null(meas_procs)) {

                # get the observation times
                obstimes <- sort(unique(unlist(lapply(meas_procs, "[[", "obstimes"))))

                # get the names of the measurement variables
                meas_vars <- unlist(lapply(meas_procs, "[[", "meas_var"))

                # initialize the matrix
                obsmat          <- matrix(0, nrow = length(obstimes), ncol = length(meas_vars)+1)
                colnames(obsmat)<- c("time", meas_vars)
                obsmat[,"time"] <- obstimes

                for(s in seq_along(meas_vars)) {
                        obsmat[, meas_vars[s]] <- ifelse(obstimes %in% meas_procs[[s]]$obstimes, 0, NA)
                }


        } else if(!is.null(datasets)) {

        }

        return(obsmat)

}