#' Return the log-likelihood of initial distribution draws
#'
#' @param initdist_objects list of initial distribution objects
#'
#' @return log likelihood
#' @export
calc_initdist_loglik = function(initdist_objects) {
    
    logliks = double(length(initdist_objects))
    
    for(s in seq_along(initdist_objects)) {
        logliks[s] = sum(dnorm(initdist_objects[[s]]$draws_cur, log = T))
    }
    
    return(sum(logliks))
}