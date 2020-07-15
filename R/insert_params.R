#' Insert natural scale parameters into a parameter matrix
#'
#' @param parmat matrix into which parameters should be inserted
#' @param param_blocks list of parameter blocks
#' @param prop should proposed parameters be inserted?
#'
#' @return insert natural scale parameters into a parameter matrix
#' @export
insert_params = function(parmat, param_blocks, prop) {
    if(prop) {
        for(s in seq_along(param_blocks)) {
            pars2parmat(parmat = parmat,
                        pars = param_blocks[[s]]$pars_prop_nat,
                        inds = param_blocks[[s]]$param_inds_nat_Cpp)
        }    
    } else {
        for(s in seq_along(param_blocks)) {
            pars2parmat(parmat = parmat,
                        pars = param_blocks[[s]]$pars_nat,
                        inds = param_blocks[[s]]$param_inds_nat_Cpp)
        }    
    }
}