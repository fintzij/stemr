#' Indicator matrix for determining which time-varying covariates change at each
#' row of a time-varying covariate matrix.
#' 
#' @param tcovar time-varying covariate matrix
#' @param tparam list for generating time-varying parameters
#' @param forcings list of forcings
#'   
#' @return indicator matrix for which time-varying covariates or parameters
#'   change at each time in a time-varying covariate matrix.
#' @export
build_tcovar_changemat <- function(tcovar, tparam = NULL, forcings = NULL) {
      
        changemat <- matrix(FALSE, nrow = nrow(tcovar), ncol = ncol(tcovar) - 1)
        colnames(changemat) <- colnames(tcovar)[-1]
        changemat[1,] <- TRUE

        for(t in 2:nrow(tcovar)) {
                changemat[t,] <- tcovar[t, 2:ncol(tcovar)] != tcovar[t-1, 2:ncol(tcovar)]
        }
        
        if(!is.null(tparam)) {
              for(s in seq_along(tparam)) {
                    changemat[,tparam[[s]]$tparam_name] <- tcovar[,1] %in% tparam[[s]]$times
              }
        }
        
        if(!is.null(forcings)) {
              for(s in seq_along(forcings)) {
                    changemat[tcovar[,forcings[[s]]$tcovar_name] == 0, forcings[[s]]$tcovar_name] <- FALSE
              }
        }

        return(changemat)
}