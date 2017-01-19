#' Construct a function for evaluating the log prior density of the initial
#' compartment counts, which are taken to have multinomial prior.
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
construct_initdist_prior <- function(state_initializer, n_strata, constants) {

        if(n_strata == 1) {
                initdist_prior_body <- paste0("dmultinom(init_counts, ",
                                              constants["popsize"],", c(",
                                              paste0(state_initializer$prior, collapse = ", "),
                                              "), log = TRUE)")

                initdist_prior <- eval(parse(text = paste0("function(init_counts) {", initdist_prior_body,"}")))

        } else {
                strata_sizes <- constants[paste0("popsize_", sapply(state_initializer,"[[","strata"))]
                initdist_prior_body <- paste0("dmultinom(init_counts[c(",
                                              apply(sapply(state_initializer,"[[","codes"), 2, paste0, collapse = ", "),
                                              ")], ",strata_sizes,", c(",
                                              apply(sapply(state_initializer,"[[","prior"), 2, paste0, collapse = ", "),
                                              "), log = TRUE)")

                initdist_prior <- eval(parse(text = paste0("function(init_counts) {sum(",
                                                           paste0(initdist_prior_body, collapse = ",\n"),")}", collapse = "\n")))
        }

        return(initdist_prior)
}