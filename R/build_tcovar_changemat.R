#' Indicator matrix for determining which time-varying covariates change at each
#' row of a time-varying covariate matrix.
#'
#' @param tcovar time-varying covariate matrix
#'
#' @return indicator matrix for which time-varying covariates change at teach
#'   time in a time-varying covariate matrix.
#' @export
build_tcovar_changemat <- function(tcovar) {
        nrow_tcovar <- nrow(tcovar); ncol_tcovar <- ncol(tcovar)
        changemat <- matrix(FALSE, nrow = nrow_tcovar, ncol = ncol_tcovar-1)
        colnames(changemat) <- colnames(tcovar)[-1]
        changemat[1,] <- TRUE

        for(t in 2:nrow(tcovar)) {
                changemat[t,] <- tcovar[t, 2:ncol_tcovar] != tcovar[t-1, 2:ncol_tcovar]
        }

        return(changemat)
}