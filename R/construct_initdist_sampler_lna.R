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
                initdist_sampler_body <- paste0("as.numeric(extraDistr::rdirichlet(1, c(",
                                              paste0(state_initializer$prior, collapse = ", "),
                                              ")))")

                initdist_sampler <- eval(parse(text = paste0("function() {", initdist_sampler_body,"}")))

        } else {
              
            # identify strata for which the initial distribution is not fixed
            not_fixed <- sapply(state_initializer, function(x) !x$fixed)
            
            # initialize the sampler body
            initdist_sampler_body <- vector("character", length(state_initializer))
            
            # construct the sampler for each stratum
            for(s in seq_along(initdist_sampler_body)) {
                  if(not_fixed[s]) {
                        
                        initdist_sampler_body[s] <- 
                              paste0("as.numeric(extraDistr::rdirichlet(1, c(",
                                     paste0(state_initializer[[s]]$prior, collapse = ", "),
                                     ")))")
                  } else {
                        initdist_sampler_body[s] <- 
                              paste0("as.numeric(c(",
                                     paste0(state_initializer[[s]]$init_states, collapse = ", "),
                                     "))")
                  }
            }
            
            initdist_sampler_body <- paste0(initdist_sampler_body, collapse = ", ")
                  
            comp_inds    <- sapply(state_initializer, "[[", "codes")
            out_order    <- order(c(comp_inds)) # order in which the counts should be returned

            initdist_sampler <- eval(parse(text = paste0("function() {c(",
                                                     paste0(initdist_sampler_body, collapse = ", "),
                                                     ")[c(",paste0(out_order, collapse = ", "),
                                                     ")]}")))
        }

        return(initdist_sampler)
}