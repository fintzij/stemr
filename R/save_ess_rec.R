#' Save the elliptical slice sampling record
#'
#' @param ess_record list with objects for recording the ESS record
#' @param ess_rec_ind C++ record index
#' @param lna_ess_schedule list with lna ESS objects
#' @param tparam list with time varying parameters
#' @param initdist_objects list with initial distribution objects
#'
#' @return inserts ESS steps and angles into the ESS record objects
#' @export
save_ess_rec = function(ess_record,
                        ess_rec_ind,
                        lna_ess_schedule = NULL,
                        tparam = NULL,
                        initdist_ess_control = NULL) {
    
    # LNA steps and angles
    if(!is.null(ess_record$lna_ess_record)) {
        mat_2_arr(dest = ess_record$lna_ess_record$ess_steps,
                  orig = cbind(sapply(lna_ess_schedule, "[[", "steps")),
                  ind = ess_rec_ind)
        
        mat_2_arr(dest = ess_record$lna_ess_record$ess_angles,
                  orig = cbind(sapply(lna_ess_schedule, "[[", "angles")),
                  ind = ess_rec_ind)
    }
    
    # initdist steps and angles
    if(!is.null(ess_record$initdist_ess_record)) {
        vec_2_mat(dest = ess_record$initdist_ess_record$ess_steps,
                  orig = initdist_ess_control$steps,
                  ind = ess_rec_ind)
        
        vec_2_mat(dest = ess_record$initdist_ess_record$ess_angles,
                  orig = initdist_ess_control$angles,
                  ind = ess_rec_ind)
    }
    
    # tparam steps and angles
    if(!is.null(ess_record$tparam_ess_record)) {
        mat_2_arr(dest = ess_record$tparam_ess_record$ess_steps,
                  orig = cbind(sapply(tparam, "[[", "steps")),
                  ind = ess_rec_ind)
        
        mat_2_arr(dest = ess_record$tparam_ess_record$ess_angles,
                  orig = cbind(sapply(tparam, "[[", "angles")),
                  ind = ess_rec_ind)
    }
    
    # increment ess_rec_ind
    increment_elem(ess_rec_ind, 0)
}
