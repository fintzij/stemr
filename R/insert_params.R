#' Insert natural scale parameters into a parameter matrix
#'
#' @param parmat matrix into which parameters should be inserted
#' @param param_blocks list of parameter blocks
#' @param nat should parameters be inserted on their estimation scale
#' @param prop should proposed parameters be inserted?
#' @param rowind C++ row index
#'
#' @return insert natural scale parameters into a parameter matrix
#' @export
insert_params = function(parmat, param_blocks, nat, prop, rowind) {
    
    if(prop) {
        if(nat) {
            for(s in seq_along(param_blocks)) {
                pars2parmat(parmat = parmat,
                            pars = param_blocks[[s]]$pars_prop_nat,
                            colinds = param_blocks[[s]]$param_inds_Cpp, 
                            rowind = rowind)
            }
        } else {
            for(s in seq_along(param_blocks)) {
                pars2parmat(parmat = parmat,
                            pars = param_blocks[[s]]$pars_prop_est,
                            colinds = param_blocks[[s]]$param_inds_Cpp, 
                            rowind = rowind)
            }
        }
    } else {
        if(nat) {
            for(s in seq_along(param_blocks)) {
                pars2parmat(parmat = parmat,
                            pars = param_blocks[[s]]$pars_nat,
                            colinds = param_blocks[[s]]$param_inds_Cpp,
                            rowind = rowind)
            }       
        } else {
            for(s in seq_along(param_blocks)) {
                pars2parmat(parmat = parmat,
                            pars = param_blocks[[s]]$pars_est,
                            colinds = param_blocks[[s]]$param_inds_Cpp,
                            rowind = rowind)
            }   
        }
    }
}
