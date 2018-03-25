#' Convert unparsed rate functions into rate functions appropriate for applying
#' the LNA to the transition event count processes.
#'
#' @param rate_fcns list of rates, as returned by stem_dynamics
#' @param compartment_codes vector of compartment codes
#' @param flow_matrix matrix indicating compartment flow, whose transpose is the
#'   stoichiometry matrix
#'
#' @return list containing a character vector with the rate functions and a
#'   named numeric vector of counting process compartment codes
#' @export
rate_fcns_4_lna <- function(rate_fcns, compartment_codes, flow_matrix) {

        # character vector for the rates
        n_rates    <- nrow(flow_matrix)
        n_comps    <- length(compartment_codes)

        # get the names of the compartments for the inflows and outflows
        if(n_rates == 1) {
                comps_from   <- rate_fcns$from
                comps_to     <- rate_fcns$to
                rate_strings <- rate_fcns$unparsed
        } else {
                comps_from   <- sapply(rate_fcns, "[[", "from")
                comps_to     <- sapply(rate_fcns, "[[", "to")
                rate_strings <- sapply(rate_fcns, "[[", "unparsed")
        }

        rate_names <- paste0(comps_from, "2", comps_to)
        comp_names <- names(compartment_codes)

        # lna_compartment codes
        lna_comp_codes        <- seq_len(n_rates) - 1
        names(lna_comp_codes) <- rate_names

        # lookup table for making substitutions
        lookup_table <- data.frame(varname     = names(compartment_codes),
                                   replacement = NA,
                                   code        = NA,
                                   stringsAsFactors = FALSE)

        lookup_table$code <- replicate(nrow(lookup_table),
                                       paste(sample(c(letters, LETTERS), 15, replace = TRUE), collapse = ""),
                                       simplify = T)

        # check that the codes do not contain the variable names
        bad_codes <- sapply(comp_names, function(x) sapply(lookup_table$code, grepl, pattern = x))
        if(any(bad_codes)) {
                for(r in seq_len(nrow(lookup_table))) {
                        while(any(sapply(comp_names, grepl, lookup_table$code[r]))) {
                                lookup_table$code[r] <- paste(sample(c(letters, LETTERS), 15, replace = TRUE), collapse = "")
                        }
                }
        }

        # convert the rates
        for(r in seq_len(n_comps)) {
              
                # get the rates that feed in and the rates that feed out
                rates_in  <- rate_names[which(flow_matrix[,r] == 1)]
                rates_out <- rate_names[which(flow_matrix[,r] == -1)]

                # generate the replacement expression
                lookup_table$replacement[r] <- paste0("(", paste0(comp_names[r],"_0"))

                if(length(rates_in)) {
                        lookup_table$replacement[r] <- paste0(lookup_table$replacement[r], " + ", paste0(rates_in, collapse = " + "))
                }

                if(length(rates_out)) {
                        lookup_table$replacement[r] <- paste0(lookup_table$replacement[r], " - ", paste0(rates_out, collapse = " - "))
                }

                lookup_table$replacement[r] <- paste0(lookup_table$replacement[r], ")")
        }

        # make the substitutions
        for(s in seq_len(n_rates)) {
                for(t in seq_len(n_comps)) {
                        rate_strings[s] <- gsub(paste0('\\<',lookup_table$varname[t],'\\>'), lookup_table$code[t], rate_strings[s])
                }
        }

        for(s in seq_len(n_rates)) {
                for(t in seq_len(n_comps)) {
                        rate_strings[s] <- gsub(paste0('\\<',lookup_table$code[t], '\\>'), lookup_table$replacement[t], rate_strings[s])
                }
        }

        return(list(lna_rates = rate_strings, lna_comp_codes = lna_comp_codes))
}