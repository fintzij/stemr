#' Construct a function for sampling new initial compartment counts, which are
#' taken to have multinomial prior.
#'
#' @param state_initializer list for intializing the comparment counts generated
#'   by \code{\link{stem_dynamics}}.
#' @param n_strata number of strata.
#' @param constants vector of model constants, which contains the population
#'   size and strata sizes.
#'
#' @return function for evaluating the log prior density of the initial
#'   compartment counts
#' @export
construct_initdist_sampler_lna <- function(state_initializer, n_strata, constants) {

        if(n_strata == 1) {
                initdist_prior_body <- paste0("as.numeric(extraDistr::rdirichlet(1, c(",
                                              paste0(state_initializer$prior, collapse = ", "),
                                              ")))")

                initdist_prior <- eval(parse(text = paste0("function() {", initdist_prior_body,"}")))

        } else {
                comp_inds    <- sapply(state_initializer, "[[", "codes")
                out_order    <- order(comp_inds) # order in which the counts should be returned
                initdist_prior_body <- paste0("as.numeric(extraDistr::rdirichlet(1, c(",
                                              apply(sapply(state_initializer,"[[","prior"), 2, paste0, collapse = ", "),
                                              ")))")

                initdist_prior <- eval(parse(text = paste0("function() {",
                                                           paste0(initdist_prior_body, collapse = ", "),"[c(",
                                                           paste0(out_order, collapse = ", "),
                                                           ")]}")))
        }

        return(initdist_prior)
}