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

                # get the observation times
                if(!is.list(datasets)) {
                        if(any(is.na(datasets))) {
                                stop("if the observation times are not the same for all of the observed variables, the datasets must be supplied as a list.")
                        }
                        obsmat <- datasets
                        colnames(obsmat) <- c("time", colnames(obsmat)[-1])

                } else {
                        # get the observation times
                        obstimes <- sort(unique(unlist(lapply(datasets, function(x) x[,1]))))

                        # get the names of the measurement variables
                        meas_vars <- sapply(lapply(datasets, function(x) x[,-1,drop=FALSE]), colnames)

                        # initialize the matrix
                        obsmat           <- matrix(NA, nrow = length(obstimes), ncol = length(meas_vars) + 1)
                        colnames(obsmat) <- c("time", meas_vars)
                        obsmat[,"time"]  <- obstimes

                        for(s in seq_along(datasets)) {
                                vars <- colnames(datasets[[s]][,-1, drop = FALSE])
                                obsmat[match(datasets[[s]][,1], obstimes), vars] <- datasets[[s]][,-1]
                        }
                }
        }

        return(obsmat)
}