#' Check if time varying parameters depend on initial conditions
#'
#' @param tparam list of time varying parameters
#' @param initdist_objects list with initial distribution objects
#' @param parmat parameter matrix
#'
#' @return time-varying parameter list with indicators for whether each
#'   parameter depends on initial conditions
#' @export
check_tpar_depends = function(tparam, initdist_objects, parmat) {
    
    # check whether initial conditions are fixed
    if(all(sapply(initdist_objects, "[[", "fixed"))) {
        
        for(s in seq_along(tparam)) {
            tparam[[s]]$update_with_inits = FALSE    
        }
        
    } else {
        # copy of parameters from matrix
        pars <- parmat[1,]
        
        # jitter initial conditions
        for(s in seq_along(initdist_objects)) {
            if(!initdist_objects[[s]]$fixed) {
                
                # new draws
                draws = rnorm(initdist_objects[[s]]$draws_cur)
                
                # calculate compartment counts
                copy_vec2(dest = pars,
                          orig = c(initdist_objects[[s]]$comp_mean +
                                       initdist_objects[[s]]$comp_sqrt_cov %*% draws),
                          inds = initdist_objects[[s]]$param_inds_Cpp)
                
                # check boundary conditions
                while(any(pars[initdist_objects[[s]]$param_inds_R] < 0) | 
                      any(pars[initdist_objects[[s]]$param_inds_R] > 
                          initdist_objects[[s]]$comp_size)) {
                    
                    # redraw and recalculate if needed
                    draws = rnorm(initdist_objects[[s]]$draws_cur)
                    
                    copy_vec2(dest = pars,
                              orig = c(initdist_objects[[s]]$comp_mean +
                                           initdist_objects[[s]]$comp_sqrt_cov %*% draws),
                              inds = initdist_objects[[s]]$param_inds_Cpp)
                }
            }
        }
        
        # now check whether the new time-varying parameters have changed
        for(s in seq_along(tparam)) {
            
            tp = tparam[[s]]$draws2par(
                    parameters = pars,
                    draws = tparam[[s]]$draws_cur)[tparam[[s]]$tpar_inds_R] 
            
            tparam[[s]]$init_dep = 
                !all.equal(tp, tparam[[s]]$tpar_cur)
        }
        
        return(tparam)
    }
}
