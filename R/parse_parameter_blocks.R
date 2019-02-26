#' Generate a list of objects used in block updating MCMC parameters
#'
#' @param parameter_blocks list of numeric or character vectors identifying
#'   parameters in each block
#' @param param_names_est character vector of (estimation scale) parameter names
#'
#' @return list containing objects used in block updates to MCMC parameters
#' @export
parse_parameter_blocks <- function(param_blocks, param_names_est) {
      
      # instatiate list
      blocks <- vector("list", length = length(param_blocks))
      
      # fill it out
      for(b in seq_along(blocks)) {
            
            # block size
            blocks[[b]]$block_size <- length(param_blocks[[b]]$pars)
            
            # indices/param names
            if(is.character(param_blocks[[b]]$pars)) {
                  
                  blocks[[b]]$param_names_est <- param_blocks[[b]]$pars
                  blocks[[b]]$param_inds_R    <- match(param_blocks[[b]]$pars, param_names_est)
                  blocks[[b]]$param_inds_Cpp  <- blocks[[b]]$param_inds_R - 1

                  if(any(is.na(blocks[[b]]$param_inds_R))) {
                        stop(paste0("Parameter name in block ", b, " does not match a parameter on the estimation scale."))
                  }
                  
            } else {
                  
                  blocks[[b]]$param_names_est <- param_names_est[param_blocks[[b]]$pars]
                  blocks[[b]]$param_inds_R    <- param_blocks[[b]]$pars
                  blocks[[b]]$param_inds_Cpp  <- blocks[[b]]$param_inds_R - 1
            }
      }
      
      # check that all parameters are accounted for
      if(!identical(sort(unlist(sapply(blocks, function(x) x$param_inds_R))), seq_along(param_names_est))) {
            stop("Not all parameters are part of a parameter block.")
      }
      
      return(blocks)
}