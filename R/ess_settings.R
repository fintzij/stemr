#' Generates a list of settings for sampling the latent LNA paths via elliptical
#' slice sampling.
#'
#' @param n_ess_updates number of elliptical slice sampling updates per MCMC
#'   iteration, defaults to 1.
#' @param ess_schedule If \code{NULL}, the LNA paths for all elementary
#'   transitions will be sampled jointly (default). May also be specified as
#'   "bystratum" if stratum specific transitions should be proposed one at a
#'   time, or as a list of charcter vectors for names of elementary transition
#'   events (i.e. each character vector corresponds to the names of rows in an
#'   LNA flow matrix) e.g., \code{"bystratum"} and \code{list(c("S_male2Imale",
#'   "I_male2R_male"), c("S_female2I_female", "I_female2R_female"))} could both
#'   be used to sample the LNA paths for males and females in a sex stratified
#'   SIR model. \code{"bystratum"} is coerced internally into the list.
#' @param randomize_order if \code{TRUE} (default) and the \code{ess_schedule}
#'   contains more than one vector, then the LNA paths for the events in the
#'   vectors are sampled in a random order in each MCMC iteration (e.g., males
#'   and females are sampled in a random order).
#' @param warmup warmup ESS updates prior to starting MCMC
#'
#' @return list with settings for elliptical slice sampling
#' @export
ess_settings <- function(n_ess_updates = 1, ess_schedule = NULL, randomize_schedule = TRUE, warmup = 200) {

        return(list(n_ess_updates = n_ess_updates,
                    ess_schedule = ess_schedule,
                    randomize_schedule = randomize_schedule,
                    warmup = warmup))

}