#' Construct a time-varying covariance matrix that includes time and
#' seasonality, based on a user-supplied time-varying covariate matrix. Used
#' internally.
#'
#' @param tcovar user-supplied time-varying covariance matrix
#' @param timestep time discretization interval
#' @param t0,tmax left and right endpoints of the time-period over which the
#'   process evolves
#'
#' @return time-varying covariate matrix that includes discretized intervals for
#'   smoothly time-varying functions. If neither a time-varying covariance
#'   matrix or a timestep is supplied, TCOVAR will just be a matrix with the
#'   left and right endpoints of the time-period.
#' @export
build_tcovar_matrix <- function(tcovar = NULL, timestep = NULL, t0, tmax, messages) {

        if(is.null(tcovar) & is.null(timestep)) {
                TCOVAR <- matrix(c(t0,tmax,t0,tmax), ncol = 2)
                colnames(TCOVAR) <- c("_time", "TIME")

        } else {
                if(!is.null(tcovar) && tcovar[1,1] != t0) {
                        warning("Time varying covariates were only specified beginning at some time after t0. It will be assumed that the values of the time-varying covariates at the first time indicated in that matrix are the same for all prior times.")
                }

                if(!is.null(timestep)) {
                        timeseq <- unique(c(seq(t0, tmax, timestep), tmax))
                } else {
                        timeseq <- NULL
                }

                TCOVAR_TIMES <- sort(unique(c(tcovar[,1, drop = FALSE], timeseq)))

                TCOVAR <- matrix(0, nrow = length(TCOVAR_TIMES), ncol = 1 + ifelse(is.null(tcovar), 0, ncol(tcovar)-1)+ ifelse(is.null(timestep), 0, 1))
                TCOVAR[,1] <- TCOVAR_TIMES

                # if there are no timevarying covariates
                if(is.null(tcovar)) {
                        TCOVAR[,2] <- TCOVAR_TIMES
                        colnames(TCOVAR) <- c("_time", "TIME")
                } else {
                        tcovar_inds <- findInterval(TCOVAR_TIMES, tcovar[,1], all.inside = T)
                        TCOVAR[,2:(ncol(TCOVAR)-1)] <- tcovar[tcovar_inds, 2:ncol(tcovar)]
                        TCOVAR[,ncol(TCOVAR)] <- timeseq[findInterval(TCOVAR_TIMES, timeseq, all.inside = TRUE)]
                        colnames(TCOVAR) <- c("_time", colnames(tcovar)[2:ncol(tcovar)], "TIME")
                }
        }

        return(TCOVAR)
}