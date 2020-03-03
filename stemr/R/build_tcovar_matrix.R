#' Construct a time-varying covariance matrix that includes time and 
#' seasonality, based on a user-supplied time-varying covariate matrix. Used 
#' internally.
#' 
#' @param tcovar user-supplied time-varying covariance matrix
#' @param tparam list for instatiating time-varying parameters
#' @param forcings list of forcings
#' @param timestep time discretization interval
#' @param census_times 
#' @param t0 
#' @param tmax 
#' @param messages 
#' @param parameters vector of model parameters
#'
#' @return matrix for time-varying covariates and parameters. If neither a
#'   time-varying covariate matrix or parameters, or a timestep is supplied,
#'   TCOVAR will just be a matrix with the left and right endpoints of the
#'   time-period.
#' @export
build_tcovar_matrix <- function(tcovar = NULL, tparam = NULL, forcings = NULL, timestep = NULL, census_times = NULL, parameters = NULL, t0, tmax, messages) {

        if(is.null(tcovar) & is.null(tparam) & is.null(timestep)) {
                TCOVAR <- matrix(c(t0,tmax,t0,tmax), ncol = 2)
                colnames(TCOVAR) <- c("_time", "TIME")

        } else {
                if(messages && !is.null(tcovar) && tcovar[1,1] != t0) {
                        warning("Time varying covariates were only specified beginning at some time after t0. It will be assumed that the values of the time-varying covariates at the first time indicated in that matrix are the same for all prior times.")
                }

                timeseq <- unique(sort(c(seq(t0, tmax, timestep), 
                                         tmax, census_times, tcovar[,1], 
                                         unlist(lapply(tparam, function(x) x$times)))))
                  
                # get times of time-varying parameters and restrict to those between t0 and tmax
                if(!is.null(tparam)) {
                      tparam_times <- sort(unique(unlist(lapply(tparam, function(x) x$times))))
                      tparam_times <- tparam_times[tparam_times >= t0 & tparam_times <= tmax]
                } else {
                      tparam_times <- NULL
                }
                
                # get tcovar times and parameters
                TCOVAR_TIMES <- sort(unique(c(tcovar[,1], tparam_times, timeseq)))

                TCOVAR <- matrix(0, nrow = length(TCOVAR_TIMES), 
                                 ncol = 1 + 
                                       ifelse(is.null(tcovar), 0, ncol(tcovar)-1) + 
                                       ifelse(is.null(tparam), 0, length(tparam)) +
                                       ifelse(is.null(timestep), 0, 1))
                TCOVAR[,1] <- TCOVAR_TIMES

                # if there are no timevarying covariates or parameters
                if(is.null(tcovar) & is.null(tparam)) {
                        TCOVAR[,2] <- TCOVAR_TIMES
                        colnames(TCOVAR) <- c("_time", "TIME")
                        
                } else {
                      
                      if(!is.null(tcovar)) {
                            tcovar_inds <- findInterval(TCOVAR_TIMES, tcovar[,1], left.open = T) + 1
                            tcovar_inds[tcovar_inds > nrow(tcovar)] <- nrow(tcovar)
                            TCOVAR[,seq(2,1 + (ncol(tcovar)-1))] <- as.matrix(tcovar[tcovar_inds, seq(2,ncol(tcovar))])
                            tcovar_names <- colnames(tcovar)[2:ncol(tcovar)]
                      } else {
                            tcovar_names <- NULL
                            tcovar_inds <- NULL
                      }
                      
                      if(!is.null(tparam)) {
                            tparam_names <- sapply(tparam, function(x) x$tparam_name)
                      } else {
                            tparam_names <- NULL
                      }
                      
                      TCOVAR[,ncol(TCOVAR)] <- timeseq[findInterval(TCOVAR_TIMES, timeseq)]
                      
                      colnames(TCOVAR) <- c("_time", tcovar_names, tparam_names, "TIME")
                }
                
                # zero out duplicated forcings
                if(!is.null(forcings)) {
                      for(f in seq_along(forcings)) {
                            forcing_times <- tcovar[which(tcovar[,forcings[[f]]$tcovar_name] != 0),1]
                            zero_inds     <- !match(round(TCOVAR_TIMES, digits = 8), round(forcing_times, digits = 8), nomatch = FALSE) 
                            TCOVAR[zero_inds,forcings[[f]]$tcovar_name] <- 0
                      }
                }
                
                # insert the time-varying parameters
                if(!is.null(tparam)) {
                      for(s in seq_along(tparam)) {
                            inds <- findInterval(TCOVAR_TIMES, tparam[[s]]$times, left.open = F)
                            inds[inds==0] <- 1
                            vals <- tparam[[s]]$draws2par(parameters = parameters, draws = tparam[[s]]$values)
                            
                            if(length(vals) != length(tparam[[s]]$times)) {
                                  stop(paste0("The draws2pars argument in tparam function number ",s," must return as many values as there are times at which the time-varying parameter is evaluated."))
                            } 
                            
                            TCOVAR[, tparam[[s]]$tparam_name] <- vals[inds]
                      }
                }
        }

        return(TCOVAR)
}