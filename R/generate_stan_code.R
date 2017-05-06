#' Generate the user defined functions for fitting the restarting,
#' log-transformed LNA for SEM counting processes via a noncentered
#' parameterization in Stan.
#'
#' @param dynamics compiled stemr dynamics
#' @param measurement_process compiled stemr measurement process
#' @param out.file optional file to output code.
#'
#' @return
generate_stan_code <- function(dynamics, measurement_process, out.file = NULL) {

        # extract objects from stem_objects
        n_times     <- measurement_process$stan_meas_objects$n_times
        n_comps     <- measurement_process$stan_meas_objects$n_comps
        n_state     <- measurement_process$stan_meas_objects$n_state
        n_emits     <- measurement_process$stan_meas_objects$n_emits
        n_modcomps  <- measurement_process$stan_meas_objects$n_modcomps
        n_params    <- measurement_process$stan_meas_objects$n_params
        hazards     <- dynamics$stan_lna_rates$hazards
        derivatives <- dynamics$stan_lna_rates$derivatives

        # derived indices
        n_modpars <- n_params - n_modcomps
        covstart  <- n_comps + 1
        covend    <- n_state
        initstart <- n_params - n_modcomps + 1
        initend   <- n_params

        # code for filling out the hazards
        HAZCODE <- paste(paste0("haz[",1:n_comps,"]=",hazards,";"), collapse = "\n")

        # code for filling out the Jacobian
        jac_inds <- matrix(c(rep(seq_len(n_comps), n_comps), rep(seq_len(n_comps), each = n_comps)), ncol = 2)
        JACCODE  <- paste(paste0("jac[",jac_inds[,1],",",jac_inds[,2],"]=",derivatives,";"), collapse = "\n")

        # read the stan function template
        tmp <- system.file(file.path("templates", "lna_stan_template.txt"), package = "stemr", mustWork = TRUE)
        con <- file(tmp)
        stan_fcns <- readLines(con)
        close(con)

        # make replacements
        stan_fcns <- gsub("__NTIMES__", n_times, stan_fcns)
        stan_fcns <- gsub("__NPARS__", n_params, stan_fcns)
        stan_fcns <- gsub("__NCOMPS__", n_comps, stan_fcns)
        stan_fcns <- gsub("__NSTATE__", n_state, stan_fcns)
        stan_fcns <- gsub("__NEMITS__", n_emits, stan_fcns)
        stan_fcns <- gsub("__NMODPARS__", n_modpars, stan_fcns)
        stan_fcns <- gsub("__NMODCOMPS__", n_modcomps, stan_fcns)
        stan_fcns <- gsub("__HAZCODE__", HAZCODE, stan_fcns)
        stan_fcns <- gsub("__JACCODE__", JACCODE, stan_fcns)
        stan_fcns <- gsub("__COVSTART__", covstart, stan_fcns)
        stan_fcns <- gsub("__COVEND__", covend, stan_fcns)
        stan_fcns <- gsub("__INITSTART__", initstart, stan_fcns)
        stan_fcns <- gsub("__INITEND__", initend, stan_fcns)

        # output to text file or return
        if(!is.null(out.file)) {
                writeLines(stan_fcns, con = out.file)
        } else {
                return(stan_fcns)
        }
}