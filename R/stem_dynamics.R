#' Generate the objects governing the dynamics of a stochastic epidemic model.
#'
#' @param rates list of rate lists, see discussion below for instructions on
#'   specifying each rate list.
#' @param parameters character vector of parameter names
#' @param compartments character vector of compartment names if there is a
#'   single stratum, or if there are multiple strata a list of character vectors
#'   where the name of each character vector is the compartment name and the
#'   character vector lists the strata to in which the compartment exists. The
#'   reserved word "ALL" can be used instead of listing all strata.
#' @param strata vector of stratum names, required only if not all compartments
#'   are common to all strata
#'
#'   (e.g. compartments = list(S = "ALL", I = "ALL", R = "ALL", D = "old");
#'   strata = c("infants", "young", "old");
#' @param constants optional character vector of constants that are referenced
#'   in the rate functions
#' @param timevar string specifying the name of the time variable, required if
#'   the time variable is referenced in one of the rates
#' @param adjacency optional matrix specifying the adjacency structure of
#'   strata, with 0 entries indicating non-adjacency and 1 for adjacency. Rows
#'   and columns must be labeled.
#'
#' @section Specifying the rate functions: Each rate function should be provided
#'   as a list with four character strings - the rate function, the compartment
#'   from which an individual exits, the compartment that she enters, and,
#'   optionally, the strata to which the rate function applies: \code{list(rate,
#'   from, to, strata = NULL)}. The \code{from} and \code{to} arguments should
#'   be strings that match entries in the \code{compartments} argument. The
#'   optional \code{strata} argument can either be a character vector of strata
#'   for which the rate applies, or may be specified as "ALL" if the rate
#'   function applies to common compartments in all strata. If there are
#'   multiple strata, the strata argument for each rate is required.
#'
#'   \emph{VERY IMPORTANT:} The rate functions must be specified at the subject
#'   level as they are parsed internally and converted to lumped rates. For
#'   example, in a standard SIR model, the subject-level infectivity and
#'   recovery rates are \eqn{\beta * I} and \eqn{\mu}, whereas the lumped rates
#'   which will be generated automatically are \eqn{\beta * I * S} and \eqn{\mu
#'   * I}.
#'
#'   The rate functions should be provided as string snippets of valid C++ code.
#'   Each function is parsed internally for references to parameters,
#'   compartments, strata, constants, and the time variable. There are a few
#'   case-sensitive reserved words that are provided to facilitate the
#'   specification of rate functions when there are multiple strata: "ALL",
#'   "ADJ", "SELF".  In an SIR model with multiple strata, for example, we could
#'   specify infectivity rates as follows: \enumerate{ \item list("beta *
#'   I_SELF", "S", "I", "ALL"): each susceptible contacts the infecteds within
#'   her own stratum at rate \eqn{\beta * I}, but does not contract an infection
#'   from outside her stratum. This rate applies to subjects in all strata.
#'   \item list("beta * pow(I_SELF, alpha)", "S", "I", "ALL"): each susceptible
#'   contacts the infecteds within her own stratum at rate \eqn{\beta *
#'   I^\alpha}, but does not contract an infection from outside her stratum.
#'   Note that the rate must be specified as valid C++, so exponentiation is
#'   done using the \code{pow()} function. \item list("beta * sum(I_ALL)", "S",
#'   "I", "ALL"): each susceptible contacts all infecteds in the population at
#'   rate \eqn{\beta * \sum_s I_s}, regardless of stratum. Note that 'I_ALL'
#'   will be replaced by a vector of I_strata, therefore using 'I_all' outside
#'   of a function, e.g. beta * I_ALL, will result in an error. However,
#'   beta*sum(I_ALL) is well defined. \item list("beta1 * I_SELF + beta2 *
#'   sum(I_adj)", "S", "I", "ALL"): each susceptible contacts infecteds in her
#'   own stratum at rate \eqn{\beta1 * I_SELF}, and contacts infecteds in
#'   adjacent strata (specified in the adjacency matrix) at rate \eqn{\beta *
#'   \sum_{s: stratum 's' is adjacent} I_s}. Note that in this case, the
#'   diagonal entries in the adjacency matrix should be set to zero.}
#'
#'   The \code{from} and \code{to} arguments will automatically be updated to
#'   reflect the strata if there are multiple strata and if no strata are
#'   specified in those arguments. Suppose, for example, that we define strata
#'   to be two adjacent regions, east and west. The first rate given above will
#'   be understood as governing transitions from S_east to I_east, and from
#'   S_west to I_west. If we we want to define rates for transitions where a
#'   subject crosses strata, we can do so by including the strata in the
#'   \code{from} and \code{to} arguments. \enumerate{\item list("beta * I",
#'   "S_east","I_west", "east"): each susceptible in the east contacts infected
#'   individuals in the east at rate \eqn{\beta * I_east}, and then transitions
#'   to the infected compartment in the west stratum. \item list("beta *
#'   I_SELF", "S", "I_west", "ALL"): all susceptibles contact infecteds in their
#'   own region at rate \eqn{\beta * I}, but upon becoming infected end up in
#'   the infected compartment in the west stratum.}
#'
#' @return list with evaluated rate functions and objects for managing the
#'   bookkeeping for epidemic paths.
#' @export
#'
stem_dynamics <- function(rates, parameters, compartments, strata = NULL, constants = NULL, timevar = NULL, adjacency = NULL) {

        # check consistency of specification and throw errors if inconsistent
        if(is.list(compartments) && ("ALL" %in% compartments) && is.null(strata)) {
                stop("In order to use 'ALL' to specify that some compartments are present for all strata, a character vector of stratum names must be supplied")
        }

        # check that none of the strata are named "ALL", "ADJ", or "SELF", and
        # that these strings do not appear in subscripts
        if(any(c("ALL", "ADJ", "SELF") %in% strata) ||
           (is.list(compartments) && any(sapply(compartments, function(x) any(c("ALL", "ADJ", "SELF") %in% x))))) {
                stop("'ALL', 'ADJ', and 'SELF' are reserved words and cannot be used in stratum or compartment names")
        }

        # check that the strata names in the compartment list match the strata
        if(is.list(compartments) &&
           any(compartments == "ALL") &&
           any(!compartments[compartments != "ALL"] %in% strata)) {
                stop(sQuote(compartments[which(!compartments %in% strata)]), "is not in the supplied vector of stratum names")
        }

        # check that the strata argument is specified for all rates if there are multiple strata
        if(!is.null(strata)) {
                for(k in 1:length(rates)) {
                        if(length(rates[[k]]) != 4) {
                                stop("If strata are specified in the model, they must also be specified in all rate functions.")
                        }
                }
        }

        # build the vector of full compartment names
        if(!is.list(compartments)) {
                compartment_names <- names(compartments)

        } else if(is.list(compartments)){

                comp_names <- names(compartments) # just the names of compartments without strata
                compartment_names <- vector("list", length(compartments)) # list with full compartment_strata names

                for(k in seq_along(compartment_names)) {
                        if(!identical(compartments[[k]], "ALL")) {
                                compartment_names[[k]] <- do.call(paste, list(comp_names[k], compartments[[k]], sep = "_"))
                        } else {
                                compartment_names[[k]] <- do.call(paste, list(comp_names[k], strata, sep = "_"))
                        }
                }
                compartment_names <- unlist(compartment_names)
        }

        # get the number of strata and the number of compartments
        n_strata <- length(strata)
        n_compartments <- length(compartment_names)

        # construct the mapping for the compartment_strata to the columns in the bookkeeping matrix
        compartment_codes <- 0:(n_compartments - 1); names(compartment_codes) <- compartment_names

        # if there are multiple strata, ensure that the rate functions are
        # specified for all compartments in their respective risk sets
        if(!is.null(strata)) {

                rate_fcns <- vector(mode = "list", length = length(rates))
                comp_all  <- paste(comp_names, "ALL", sep = "_")
                comp_adj  <- paste(comp_names, "ADJ", sep = "_")

                for(k in seq_along(rates)) {

                        # get the relevant strata
                        if(identical(rates[[k]][[4]], "ALL")) {
                                rel_strata <- strata
                        } else {
                                rel_strata <-rates[[k]][[4]]
                        }

                        # ensure that there is one rate function per stratum
                        rate_fcns[[k]] <- vector(mode = "list", length = length(rel_strata))

                        # modify the rate function to reflect the stratum
                        for(j in seq_along(rel_strata)) {
                                strat_self <- rel_strata[j] # get stratum name

                                rate_fcns[[k]][[j]] <- vector(mode = "list", length = 4)
                                rate_fcns[[k]][[j]][c(1,3,4)] <- rates[[k]][1:3] # copy generic rate function

                                # make substitutions in the rate string
                                # SELF substitution
                                if(grepl("SELF", rates[[k]][[1]])) {
                                        rate_fcns[[k]][[j]][[1]] <- gsub("SELF", strat_self, rate_fcns[[k]][[j]][[1]])
                                }

                                # ALL substitution
                                if(grepl("ALL", rates[[k]][[1]])) {
                                        which_all <- sapply(comp_all, FUN = grepl, rate_fcns[[k]][[j]][[1]])
                                        which_sub <- comp_all[which_all]

                                        for(l in seq_along(which_sub)) {
                                                sub_comp <- comp_names[which(comp_all == which_sub[l])]
                                                sub_all  <- compartment_names[compartment_names %in% paste(sub_comp, strata, sep = "_")]
                                                rate_fcns[[k]][[j]][[1]] <- gsub(which_sub[l], paste(sub_all, collapse = ", "), rate_fcns[[k]][[j]][[1]])
                                        }
                                }

                                # ADJ substitution
                                if(grepl("ADJ", rates[[k]][[1]])) {
                                        # vector of adjacent strata
                                        strat_adj <- colnames(adjacency)[adjacency[strat_self,] == 1]

                                        which_adj <- sapply(comp_adj, FUN = grepl, rate_fcns[[k]][[j]][[1]])
                                        which_sub <- comp_adj[which_adj]

                                        for(l in seq_along(which_sub)) {
                                                sub_comp <- comp_names[which(comp_adj == which_sub[l])]
                                                sub_adj  <- compartment_names[grepl(pattern = sub_comp, x = compartment_names)]
                                                sub_adj  <- compartment_names[which(compartment_names %in% paste(sub_comp, strat_adj, sep = "_"))]

                                                rate_fcns[[k]][[j]][[1]] <- gsub(which_sub[l], paste(sub_adj, collapse = ", "), rate_fcns[[k]][[j]][[1]])
                                        }
                                }

                                # make substitutions in the 'from' and 'to' arguments
                                if((!rate_fcns[[k]][[j]][[3]] %in% compartment_names) && (rate_fcns[[k]][[j]][[3]] %in% comp_names)) {
                                        rate_fcns[[k]][[j]][[3]] <- paste(rate_fcns[[k]][[j]][[3]], strat_self, sep = "_")
                                }

                                if((!rate_fcns[[k]][[j]][[4]] %in% compartment_names) && (rate_fcns[[k]][[j]][[4]] %in% comp_names)) {
                                        rate_fcns[[k]][[j]][[4]] <- paste(rate_fcns[[k]][[j]][[4]], strat_self, sep = "_")
                                }

                                # instatiate lumped rate functions
                                rate_fcns[[k]][[j]][[2]] <- paste0("(", rate_fcns[[k]][[j]][[1]], ") * ", rate_fcns[[k]][[j]][[3]])
                        }
                }

                rate_fcns <- unlist(rate_fcns, recursive = FALSE)

        } else {
                rate_fcns <- rates
        }

        # construct the flow matrix
        flow_matrix <- build_flowmat(rate_fcns, compartment_names)

        # parse the rates ######### LEFT OFF HERE. FINISHED FUNCTION TO BUILD FLOW MATRIX. NEXT, FUNCTION TO BUILD THE RATE ADJACENCY MATRIX AND THE C++ FUNCTIONS TO COMPILE AND CALL THE RATE FUNCTIONS ######
        rate_vec <- vector(mode = "list", )
        for(k in seq_along(rate_fcns)) {
                rate_vec[k] <- parse_rate(rates[[k]], parameters, compartment_names, constants, strata, adjacency)
        }

}