#' Reformat the rates for use in the LNA and compute the derivatives of the rate functions.
#'
#' @param rates list of rates, as returned by stem_dynamics
#' @param parameters vector of parameter codes
#' @param tcovar vector of time-varying covariate codes
#' @param constants vector of constant codes
#' @param compartments vector of compartment codes
#'
#' @return character vector of rate derivatives ordered by compartment within
#'   rates
#' @export
get_lna_rates <- function(rates, parameters, tcovar, constants, compartments) {

        # Underscores will need to be removed and re-inserted at the end.
        # Thus, we create a lookup table for the strings with and without underscores
        lookup_table      <- data.frame(varname     = c(paste("parameters[",0:(length(parameters)-1), "]", sep = ""),
                                                        paste("constants[", 0:(length(constants)-1), "]", sep = ""),
                                                        paste("state[",0:(length(compartments)-1), "]", sep = ""),
                                                        paste("tcovar[", 1:(length(tcovar)), "]", sep = "")),
                                        search_name = c(paste("parameters\\[",0:(length(parameters)-1), "\\]", sep = ""),
                                                        paste("constants\\[",0:(length(constants)-1), "\\]", sep = ""),
                                                        paste("state\\[",0:(length(compartments)-1), "\\]", sep = ""),
                                                        paste("tcovar\\[",1:(length(tcovar)), "\\]", sep = "")),
                                        code        = NA,
                                        stringsAsFactors = FALSE)
        lookup_table$code <- replicate(nrow(lookup_table), paste(sample(c(letters, LETTERS), 15, replace = TRUE), collapse = ""), simplify = T)

        if("TIME" %in% names(tcovar)) lookup_table[which(lookup_table[,2] == paste("tcovar\\[",which(names(tcovar) == "TIME"),"\\]", sep = "")), 1] <- "t"

        # get indices for which rows correspond to the compartments
        comp_inds <- unname(sapply(paste("state\\[",0:(length(compartments)-1), "\\]", sep = ""), match, table = lookup_table[,2]))

        # get the unparsed rate strings
        rate_strings <- sapply(rates, "[[", "lumped")

        # make the substitutions from the lookup table
        for(s in seq_along(rate_strings)) {
                for(j in seq_len(nrow(lookup_table))) {
                        rate_strings[s] <- gsub(pattern = lookup_table[j,2], replacement = lookup_table[j,3], x = rate_strings[s])
                }
        }

        # generate symbolic expressions for the rates and other objects
        rate_syms <- lapply(rate_strings, function(x) (parse(text = x)))
        comp_syms <- lapply(lookup_table[comp_inds,3], Ryacas::Sym)

        # create list to store derivatives
        derivatives <- vector(mode = "list", length = length(rate_strings))

        # compute the derivatives
        for(t in seq_along(rate_strings)) {
                derivatives[[t]] <- vector(mode = "list", length = length(compartments))

                for(s in seq_along(derivatives[[t]])) {
                        derivatives[[t]][[s]] <- D(rate_syms[[t]], comp_syms[[s]])
                        derivatives[[t]][[s]] <- paste(deparse(derivatives[[t]][[s]]), collapse = "")
                }

        }
        derivatives   <- unlist(derivatives)

        # replace the hash codes with the names of the vector elements
        for(s in seq_along(rate_strings)) {
                for(j in seq_len(nrow(lookup_table))) {
                        rate_strings[s] <- gsub(pattern = lookup_table[j,3], replacement = lookup_table[j,1], x = rate_strings[s])
                        rate_strings[s] <- gsub(" ", "", rate_strings[s])
                }
        }

        for(s in seq_along(derivatives)) {
                for(j in seq_len(nrow(lookup_table))) {
                        derivatives[s] <- gsub(pattern = lookup_table[j,3], replacement = lookup_table[j,1], x = derivatives[s])
                        derivatives[s] <- gsub(" ", "", derivatives[s])
                }
        }

        return(list(hazards = rate_strings, derivatives = derivatives))
}