#' Insert natural scale parameters into a parameter matrix
#'
#' @param parmat matrix into which parameters should be inserted
#' @param initdist_objects list of parameter blocks
#' @param prop should proposed parameters be inserted?
#' @param rowind row index
#'
#' @return insert natural scale parameters into a parameter matrix
#' @export
insert_initdist = function(parmat, initdist_objects, prop = FALSE, rowind = 0) {
    if(prop) {
        for(s in seq_along(initdist_objects)) {
            pars2parmat(parmat = parmat,
                        pars = initdist_objects[[s]]$init_volumes_prop,
                        colinds = initdist_objects[[s]]$param_inds_Cpp,
                        rowind = rowind)
        }    
    } else {
        for(s in seq_along(initdist_objects)) {
            pars2parmat(parmat = parmat,
                        pars = initdist_objects[[s]]$init_volumes,
                        colinds = initdist_objects[[s]]$param_inds_Cpp,
                        rowind = rowind)
        }    
    }
}