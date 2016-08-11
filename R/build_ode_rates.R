#' Construct the ODE hazard functions
#'
#' @inheritParams build_ode_rates
#'
#' @return list of rate function strings and parameter codes
#' @export
build_ode_rates <- function(rate_fcns, param_codes, const_codes, tcovar_codes, compartment_codes) {

        ode_param_codes <- c(param_codes, const_codes + length(param_codes), tcovar_codes + length(param_codes) + length(const_codes)-1)

        # Underscores will need to be removed and re-inserted at the end.
        # Thus, we create a lookup table for the strings with and without underscores
        lookup_table      <- data.frame(varname     = c(paste("pars[", ode_param_codes, "]", sep = ""),
                                                        paste("x[", compartment_codes, "]", sep = "")),
                                        search_name = c(paste("parameters\\[", param_codes, "\\]", sep = ""),
                                                        paste("constants\\[", const_codes, "\\]", sep = ""),
                                                        paste("tcovar\\[", tcovar_codes, "\\]", sep = ""),
                                                        paste("state\\[", compartment_codes, "\\]", sep = "")),
                                        stringsAsFactors = FALSE)

        if("TIME" %in% names(tcovar_codes)) lookup_table[which(lookup_table[,2] == paste("tcovar\\[",which(names(tcovar_codes) == "TIME"),"\\]", sep = "")), 1] <- "t"

        # get indices for which rows correspond to the compartments
        comp_inds <- unname(sapply(paste("state\\[",0:(length(compartments)-1), "\\]", sep = ""), match, table = lookup_table[,2]))

        # get the unparsed rate strings
        rate_strings <- sapply(rate_fcns, "[[", "lumped")

        # make the substitutions from the lookup table
        for(s in seq_along(rate_strings)) {
                for(j in seq_len(nrow(lookup_table))) {
                        rate_strings[s] <- gsub(pattern = lookup_table[j,2], replacement = lookup_table[j,1], x = rate_strings[s])
                }
        }

        return(list(hazards = rate_strings, ode_param_codes = ode_param_codes))
}