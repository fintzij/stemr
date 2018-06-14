#' Reconstitute a joint kernel covariance matrix from block covariances
#'
#' @param kernel_objects list of kernel objects containing the sub-block
#'   covariance matrices
#' @param parameter_blocks list of parameter blocks
#'
#' @return reconstituted covariance matrix
#' @export
blocks2cov <- function(kernel_objects, parameter_blocks) {
      
      n_params <- sum(sapply(parameter_blocks, function(x) x$block_size))
      covmat   <- matrix(0.0, n_params, n_params)
      param_names <- vector("character", n_params)
      
      # insert blocks
      for(b in seq_along(parameter_blocks)) {
            covmat[parameter_blocks[[b]]$param_inds_R, parameter_blocks[[b]]$param_inds_R] <- 
                  kernel_objects[[b]]$kernel_cov
            param_names[parameter_blocks[[b]]$param_inds_R] <- parameter_blocks[[b]]$param_names_est
      }
      
      # set row and column names
      rownames(covmat) <- colnames(covmat) <- param_names
      
      # return covmat
      return(covmat)
}