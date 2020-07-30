#' Record ESS steps and angles
#'
#' @param ess_record object where the steps and angles are to be stored
#' @param ess_rec_ind index
#' @param lna_ess_schedule lna ESS list
#' @param tparam tparam ESS list
#' @param initdist_objects 
#'
#' @return save the ESS record
save_ess_record = 
    function(ess_record, 
             ess_rec_ind,
             lna_ess_schedule = NULL,
             initdist_objects = NULL,
             tparam = NULL) {
        
        # LNA ESS record
        if(!is.null(ess_record$lna_ess_record)) {
            for(s in seq_along(lna_ess_schedule)) {
                vec_2_arr(dest = ess_record$lna_ess_record$ess_steps,
                          orig = lna_ess_schedule[[s]]$steps,
                          col_ind = s-1,
                          slice_ind = ess_rec_ind)
                
                vec_2_arr(dest = ess_record$lna_ess_record$ess_angles,
                          orig = lna_ess_schedule[[s]]$angles,
                          col_ind = s-1,
                          slice_ind = ess_rec_ind)
            }
        }
        
        # initdist ESS record
        if(!is.null(ess_record$initdist_ess_record)) {
            vec_2_mat(dest = ess_record$initdist_ess_record$ess_steps,
                      orig = initdist_objects$initdist_steps,
                      ind = ess_rec_ind)
            
            vec_2_mat(dest = ess_record$initdist_ess_record$ess_angles,
                      orig = initdist_objects$initdist_angle,
                      ind = ess_rec_ind)
        }
        
        # tparam ESS record
        if(!is.null(ess_record$tparam_ess_record)) {
            for(s in seq_along(tparam)) {
                vec_2_arr(dest = ess_record$tparam$ess_steps,
                          orig = tparam[[s]]$steps,
                          col_ind = s-1,
                          slice_ind = ess_rec_ind)
                
                vec_2_arr(dest = ess_record$tparam$ess_angles,
                          orig = tparam[[s]]$angles,
                          col_ind = s-1,
                          slice_ind = ess_rec_ind)
            }
        }
        
        # increment rec_ind
        increment_vec(ess_rec_ind, 1)
    }
