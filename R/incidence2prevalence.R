#' Convert an LNA/ODE incidence path to a prevalence path.
#'
#' @param path matrix containing the incidence path the first column of which
#'   are incidence times
#' @param flow_matrix flow matrix
#' @param init_state vector of initial compartment counts
#' @param forcings list of forcing lists
#' @param forcing_matrix matrix containing the forcings the first column of
#'   which are times of forcings
#'
#' @return matrix containing the prevalence path
#' @export
incidence2prevalence <- function(path, flow_matrix, init_state, forcings = NULL, forcing_matrix = NULL) {
      
      # make sure that the objects are of the appropriate type
      if(class(path) != "matrix") path <- as.matrix(path)
      if(class(flow_matrix) != "matrix") flow_matrix <- as.matrix(flow_matrix)
      if(class(init_state) != "numeric") init_state  <- as.numeric(init_state)
      if(!is.null(forcing_matrix) && class(forcing_matrix) != "matrix") forcing_matrix <- as.matrix(forcing_matrix)
      
      if(!is.null(forcings)) {
            
            # get times
            path_times    <- path[,1]
            forcing_times <- forcing_matrix[,1]
            
            if(length(path_times) != length(forcing_times)) {
                  
                  # all times, sorted
                  census_times <- sort(unique(c(path_times, forcing_times)))
                  
                  # indices of path and forcing indices
                  path_inds    <- census_times %in% path_times
                  forcing_inds <- census_times %in% forcing_times
                  
                  # create matrices and insert times
                  path_mtx     <- matrix(0.0, nrow = length(census_times), ncol = ncol(path))
                  forcing_mtx  <- matrix(0.0, nrow = length(census_times), ncol = ncol(forcing_matrix))
                  path_mtx[,1] <- census_times
                  forcing_mtx[,1] <- census_times
                  
                  # insert path and forcings
                  path_mtx[path_inds,-1]       <- path[,-1]
                  forcing_mtx[forcing_inds,-1] <- forcing_matrix[,-1]
                  
            } else {
                  path_mtx     <- path
                  forcing_mtx  <- forcing_matrix
                  forcing_inds <- rep(TRUE, nrow(path))
            }
            
            colnames(path_mtx)   <- colnames(path)
            colnames(forcing_mtx) <- colnames(forcing_matrix)
            
            # names and indices
            forcing_tcovars   <- sapply(forcings, function(x) x$tcovar_name)
            forcing_tcov_inds <- match(forcing_tcovars, colnames(forcing_matrix)) - 1
            forcing_events    <- c(sapply(forcings, function(x) paste0(x$from, "2", x$to)))
            
            # matrix indicating which compartments are involved in which forcings in and out
            forcings_out <- matrix(0.0, 
                                   nrow = ncol(flow_matrix), ncol = length(forcings),
                                   dimnames = list(colnames(flow_matrix), forcing_tcovars))
            
            forcing_transfers <- array(0.0, 
                                       dim = c(ncol(flow_matrix),
                                               ncol(flow_matrix),
                                               length(forcings)),
                                       dimnames = list(colnames(flow_matrix),
                                                       colnames(flow_matrix),
                                                       forcing_tcovars))
            
            for(s in seq_along(forcings)) {
                  
                  forcings_out[forcings[[s]]$from, s] <- 1
                  
                  for(t in seq_along(forcings[[s]]$from)) {
                        forcing_transfers[forcings[[s]]$from[t], forcings[[s]]$from[t], s] <- -1
                        forcing_transfers[forcings[[s]]$to[t], forcings[[s]]$from[t], s]    <- 1
                  }
            }
            
      } else {
            path_mtx          <- path
            forcing_inds      <- rep(FALSE, nrow(path))
            forcing_tcovars   <- character(0L)
            forcing_tcov_inds <- integer(0L)
            forcing_events    <- character(0L)
            forcings_out      <- matrix(0.0, nrow = 0, ncol = 0)
            forcing_transfers <- array(0.0, dim = c(0,0,0))
      }
      
      conv_path <-
            lna_incid2prev(
                  path              = path_mtx,
                  flow_matrix       = flow_matrix,
                  init_state        = init_state,
                  forcing_matrix    = forcing_mtx,
                  forcing_inds      = forcing_inds,
                  forcing_tcov_inds = forcing_tcov_inds,
                  forcings_out      = forcings_out,
                  forcing_transfers = forcing_transfers
            )
      
      if(!is.null(forcings)) {
            conv_path <- conv_path[conv_path[,1] == path_times,]
      }
      
      timename <- ifelse(is.null(colnames(path)), "time", colnames(path)[1]) 
      colnames(conv_path) <- c(timename, colnames(flow_matrix))
      
      return(conv_path)
}