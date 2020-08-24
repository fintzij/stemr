#' Prepare an LNA elliptical slice sampling schedule
#'
#' @param stem_object stem_object list
#' @param initializer initializer list
#' @param lna_ess_control LNA control list
#'
#' @return LNA ESS schedule
#' @export
prepare_lna_ess_schedule = function(stem_object, initializer, lna_ess_control) {
    
    ess_schedule = 
        if(lna_ess_control$joint_strata_update) {
            vector("list", length = 1)
        } else {
            vector("list", length = stem_object$dynamics$n_strata)
        }
    
    # indices
    if(length(ess_schedule) == 1) {
        
        # all inds and initdist codes
        ess_schedule[[1]]$ess_inds = seq_len(nrow(stem_object$dynamics$flow_matrix))
        ess_schedule[[1]]$initdist_codes = seq_along(initializer)
        
        # complementary inds
        ess_schedule[[1]]$complementary_inds <- 
            setdiff(seq_len(nrow(stem_object$dynamics$flow_matrix)), 
                    ess_schedule[[1]]$ess_inds)
        
        # stratum code
        ess_schedule[[1]]$stratum_code = 1
        
        # control
        ess_schedule[[1]]$bracket_width = lna_ess_control$bracket_width
        ess_schedule[[1]]$bracket_update_iter = lna_ess_control$bracket_update_iter
        ess_schedule[[1]]$angle_mean  = 0
        ess_schedule[[1]]$angle_var   = pi^2/3
        ess_schedule[[1]]$angle_resid = 0
        ess_schedule[[1]]$n_updates   = lna_ess_control$n_updates 
        
    } else {
        ess_inds =
            lapply(paste0("_", names(stem_object$dynamics$strata_codes)),
                   function(x) grep(x, rownames(stem_object$dynamics$flow_matrix)))
        
        initdist_codes = 
            lapply(names(ess_schedule$strata_codes),
                   function(x) match(x, sapply(initializer, "[[", "strata")))
        
        complementary_inds <- 
            lapply(ess_inds, 
                   function(x)
                       setdiff(seq_len(nrow(stem_object$dynamics$flow_matrix)), x))
        
        for(s in seq_along(ess_schedule)) {
            ess_schedule[[s]]$ess_inds = ess_inds[[s]]
            ess_schedule[[s]]$initdist_codes = initdist_codes[[s]]
            ess_schedule[[s]]$complementary_inds = complementary_inds[[s]]
            ess_schedule[[s]]$stratum_code = s
            ess_schedule[[s]]$bracket_width = lna_ess_control$bracket_width
            ess_schedule[[s]]$bracket_update_iter = lna_ess_control$bracket_update_iter
            ess_schedule[[s]]$angle_mean  = 0
            ess_schedule[[s]]$angle_var   = pi^2/3
            ess_schedule[[s]]$angle_resid = 0
            ess_schedule[[s]]$n_updates   = lna_ess_control$n_updates 
        }
    }
    
    # return the ess_schedule object
    return(ess_schedule)
}
