#' Generate the user defined functions for fitting the restarting,
#' log-transformed LNA for SEM counting processes via a noncentered
#' parameterization in Stan.
#'
#' @param dynamics compiled stemr dynamics
#' @param measurement_process compiled stemr measurement process
#' @param separate_params should the model parameters be separated out (default) in the
#'   template or combined (e.g. for multivariate normal) for convenience?
#' @param out.file optional file to output code.
#'
#' @return
generate_stan_code <- function(dynamics, measurement_process, separate_params = TRUE, out.file = NULL) {

        # extract objects from stem_objects
        n_times     <- measurement_process$stan_meas_objects$n_times
        n_comps     <- measurement_process$stan_meas_objects$n_comps
        n_state     <- measurement_process$stan_meas_objects$n_state
        n_emits     <- measurement_process$stan_meas_objects$n_emits
        n_modcomps  <- measurement_process$stan_meas_objects$n_modcomps
        n_params    <- measurement_process$stan_meas_objects$n_params
        emitinds    <- measurement_process$obscomp_codes
        hazards     <- dynamics$stan_lna_rates$hazards
        derivatives <- dynamics$stan_lna_rates$derivatives
        initstates  <- paste(colnames(dynamics$flow_matrix_lna), collapse = "\n")
        lna_param_codes <- dynamics$stan_lna_rates$lna_param_codes

        # derived indices
        n_modpars <- n_params - n_modcomps
        covstart  <- n_comps + 1
        covend    <- n_state
        initstart <- n_params - n_modcomps + 1
        initend   <- n_params

        # code for filling out the hazards
        HAZCODE <- paste(paste0("haz[",1:n_comps,"]=",hazards,";"), collapse = "\n")

        # code for filling out the Jacobian
        nonzeros <- which(derivatives != "0")
        jac_inds <- matrix(c(rep(seq_len(n_comps), n_comps), rep(seq_len(n_comps), each = n_comps)), ncol = 2)
        JACCODE  <- paste(paste0("jac[",jac_inds[nonzeros,1],",",jac_inds[nonzeros,2],"]=",
                                 derivatives[nonzeros],";"), collapse = "\n")

        # code for inserting a template for the raw model parameters, transformed parameters, and priors
        if(separate_params) {
                RAWPARAMS <- paste(paste0("real ", names(lna_param_codes)[seq_len(n_modpars)],"_raw;"), collapse = "\n")
                NATPARAMS <- paste(paste0("real ", names(lna_param_codes)[seq_len(n_modpars)],";"), collapse = "\n")
                PARAMPRIORS <- paste(paste0(names(lna_param_codes)[seq_len(n_modpars)],"_raw ~ __PRIORDIST__;"), collapse = "\n")
        } else {
                RAWPARAMS <- paste(paste0("// theta_raw[",seq_len(n_modpars),"] = ",
                                          names(lna_param_codes)[seq_len(n_modpars)],"_raw;"), collapse = "\n")
                RAWPARAMS <- paste(RAWPARAMS, paste0("vector[",n_modpars,"] theta_raw;"), sep = "\n")

                NATPARAMS <- paste(paste0("// theta[",seq_len(n_modpars),"] = ",
                                          names(lna_param_codes)[seq_len(n_modpars)],";"), collapse = "\n")
                NATPARAMS <- paste(NATPARAMS, paste0("real theta[",n_modpars,"];"), sep = "\n")

                PARAMPRIORS <- "theta_raw ~ __PRIORDIST__;"
        }

        # population size constants
        which_popsizes <- grep("popsize", names(dynamics$constants))
        popsizes       <- paste(paste0("real ", names(dynamics$constants)[which_popsizes]," = ",
                                        dynamics$constants[which_popsizes], ";"), collapse = "\n")
        log_popsizes   <- paste(paste0("real log_", names(dynamics$constants)[which_popsizes]," = log(",
                                        names(dynamics$constants)[which_popsizes], ");"), collapse = "\n")
        popsizes       <- paste(popsizes, log_popsizes, sep = "\n")

        # strata indices in the initial state vector
        if(dynamics$n_strata == 1) {
                strata_inds <- ""
        } else {
                strata_inds <- vector("character", length = dynamics$n_strata)
                for(k in seq_along(strata_inds)) {
                        strata_inds[k] <- paste0("int ",dynamics$state_initializer[[k]]$strata,
                                                 "_inds[",length(dynamics$state_initializer[[k]]$codes),"] = {",
                                                 paste0(dynamics$state_initializer[[k]]$codes + n_modpars, collapse=","),"};")
                }
                strata_inds <- paste(strata_inds, collapse = "\n")
        }

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
        stan_fcns <- gsub("__RAWPARAMS__", RAWPARAMS, stan_fcns)
        stan_fcns <- gsub("__NATPARAMS__", NATPARAMS, stan_fcns)
        stan_fcns <- gsub("__PARAMPRIORS__", PARAMPRIORS, stan_fcns)
        stan_fcns <- gsub("__EMITINDS__", paste(emitinds, collapse = ","), stan_fcns)
        stan_fcns <- gsub("__INITSTATES__", initstates, stan_fcns)
        stan_fcns <- gsub("__POPSIZES__", popsizes, stan_fcns)
        stan_fcns <- gsub("__STRATAINDS__", strata_inds, stan_fcns)

        # output to text file or return
        if(!is.null(out.file)) {
                writeLines(stan_fcns, con = out.file)
        } else {
                return(stan_fcns)
        }
}