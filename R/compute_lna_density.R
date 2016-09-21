#' Compute the log-transition density for a latent LNA path
#'
#' @param path list containing the current path and the paths of its ODEs
#' @param fixed_inits is the initial state vector fixed?
#' @param log_scale is the LNA applied to the log transformed process?
#' @param incidence_codes vector of incidence codes
#'
#' @return log density of the sampled LNA path
#' @export
compute_lna_density <- function(path, fixed_inits, log_scale, incidence_codes) {

        n_comps <- ncol(path$drift_process)
        n_times <- nrow(path$drift_process)

        if(!is.null(incidence_codes)) {
                incid_codes_step <- incidence_codes + 1
                incid_codes_path <- incidence_codes + 2
                n_prev  <- n_comps - length(incidence_codes)
        }

        if(fixed_inits) {
                .log_lik <- 0
        } else {
                .log_lik <- tmvtnorm::dtmvnorm(x = path$path[1,-1],
                                               mean = path$drift_process[1,],
                                               sigma = path$diffusion_process[,,1],
                                               lower = ifelse(log_scale, rep(-Inf, n_comps), rep(0, n_comps)),
                                               log = T)
        }

        if(log_scale) {
                for(j in 2:n_times) {
                        if(!is.null(incidence_codes)) {
                                .log_lik <- .log_lik + tmvtnorm::dtmvnorm(path$path[j,-1],
                                                                          path$drift_process[j,] + path$residual_process[j,],
                                                                          path$diffusion_process[,,j],
                                                                          lower = c(rep(-Inf, n_prev), path$path[j-1, incid_codes_path]),
                                                                          log = TRUE)
                        } else {
                                .log_lik <- .log_lik + tmvtnorm::dtmvnorm(path$path[j,-1],
                                                                          path$drift_process[j,] + path$residual_process[j,],
                                                                          path$diffusion_process[,,j], log = TRUE)
                        }
                }

        } else if(!log_scale){
                for(j in 2:n_times) {
                        if(!is.null(incidence_codes)) {
                                .log_lik <- .log_lik + tmvtnorm::dtmvnorm(path$path[j,-1],
                                                                          path$drift_process[j,] + path$residual_process[j,],
                                                                          path$diffusion_process[,,j],
                                                                          lower = c(rep(0, n_prev), path$path[j-1, incid_codes_path]),
                                                                          log = TRUE)
                        } else {
                                .log_lik <- .log_lik + tmvtnorm::dtmvnorm(path$path[j,-1],
                                                                          path$drift_process[j,] + path$residual_process[j,],
                                                                          path$diffusion_process[,,j],
                                                                          lower = rep(0, n_comps), log = TRUE)
                        }
                }
        }

        return(.log_lik)
}