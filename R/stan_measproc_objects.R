#' Generates additional objects for use in constructing LNA user defined functions for use in Stan.
#'
#' @param dynamics list of parsed dynamics
#' @param meas_procs list of parsed measurment process functions
#' @param obstimes list of vectors of observation times
#'
#' @return list of additional objects for Stan user defined functions
#' @export
stan_measproc_objects <- function(dynamics, meas_procs, obsmat) {

        n_times    = nrow(obsmat)      # number of observation (census) times
        n_emits    = ncol(obsmat) - 1  # number of transition events for which an emission process is defined
        n_comps    = nrow(dynamics$flow_matrix_lna) # number of LNA compartments = # event types
        n_modcomps = ncol(dynamics$flow_matrix_lna) # number of model compartments
        n_state    = n_comps + n_comps^2 # length of the LNA state vector
        n_params   = length(dynamics$stan_lna_rates$lna_param_codes) # parameters (including initdist pars)
        n_modpars  = n_params - n_modcomps
        emit_inds  = match(names(meas_procs$incidence_codes_lna), rownames(dynamics$flow_matrix_lna))

        return(list(n_times    = n_times,
                    n_emits    = n_emits,
                    n_comps    = n_comps,
                    n_modcomps = n_modcomps,
                    n_state    = n_state,
                    n_params   = n_params,
                    n_modpars  = n_modpars,
                    emit_inds  = emit_inds))
}