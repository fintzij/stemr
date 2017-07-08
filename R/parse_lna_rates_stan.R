#' Parse the LNA rates for use in generating Stan code.
#'
#' @param lna_rates character vector of LNA rates
#' @param param_codes named numeric vector of parameter codes
#' @param const_codes named numeric vector of constant codes
#' @param tcovar_codes named numeric vector of time-varying covariate codes
#' @param lna_comp_codes named numeric vector of LNA compartment codes
#'
#' @return string snippets for the LNA
#' @export
parse_lna_rates_stan <- function(lna_rates, param_codes, const_codes, tcovar_codes, lna_comp_codes) {

        # adjust indexing for R/Stan style indexing
        param_codes  <- param_codes + 1
        const_codes  <- const_codes + 1

        # identify which constants are initial compartment counts
        which_const_inits <- grep("_0", names(const_codes))

        # concatenate the parameter codes and the initial compartment codes
        lna_param_codes <- c(param_codes, length(param_codes) + seq_along(const_codes[which_const_inits]))

        # generate lookup table
        lookup_table <- data.frame(varname     = c(paste("lna_params[", lna_param_codes, "]", sep = ""),
                                                   paste("Z[", lna_comp_codes, "]", sep = ""),
                                                   "t"),
                                   search_name = c(names(param_codes),
                                                   names(const_codes)[which_const_inits],
                                                   names(lna_comp_codes),
                                                   names(tcovar_codes)),
                                   code        = NA,
                                   log_code    = NA,
                                   stringsAsFactors = FALSE)

        # get indices for which rows correspond to the compartments
        lookup_table[which(lookup_table[,"search_name"] == "TIME"), 1] <- "t"

        comp_inds <- unname(sapply(names(lna_comp_codes), match, table = lookup_table[,"search_name"]))
        time_ind  <- match("t", lookup_table[,1])

        # generate the code strings
        lookup_table$code <- lookup_table$log_code <- replicate(nrow(lookup_table),
                                       paste(sample(c(letters, LETTERS), 15, replace = TRUE), collapse = ""),
                                       simplify = T)
        lookup_table$log_code[comp_inds] <- paste0("exp(", lookup_table$log_code[comp_inds],")")

        # grab the rates
        rates <- lna_rates

        # make the substitutions in the rate strings
        for(s in seq_along(rates)) {
                for(j in seq_len(nrow(lookup_table))) {
                        rates[s] <- gsub(pattern = lookup_table[j,"search_name"],
                                             replacement = lookup_table[j,"log_code"], x = rates[s])
                }

                for(j in seq_along(comp_inds)) {
                        rates[s] <- gsub(pattern = lookup_table[comp_inds[j], "code"],
                                             replacement = lookup_table[comp_inds[j], "log_code"], x = rates[s])
                }
        }

        # generate hazards and derivatives for the Jacobian
        hazards     <- vector(mode = "character", length = length(lna_comp_codes))
        derivatives <- vector(mode = "list", length = length(hazards))

        for(r in seq_along(hazards)) {
                hazards[r] <- paste0("(exp(-",
                                     lookup_table[comp_inds[r], "code"],
                                     ") - 0.5*exp(-2*",
                                     lookup_table[comp_inds[r], "code"],
                                     "))*(",rates[r],")")

                hazards[r] <- gsub("\\+ \\-", "- ", hazards[r])
        }

        # generate symbolic expressions for the rates and other objects
        rate_syms   <- lapply(hazards, function(x) (parse(text = x)))
        comp_syms   <- lapply(lookup_table[comp_inds,"code"], Ryacas::Sym)
        time_derivs <- vector(mode = "list", length = length(hazards))

        # if there are time-varying rates, generate the time-derivative symbols
        if(!is.na(time_ind)) time_sym  <- Ryacas::Sym(lookup_table[time_ind,"code"])

        # compute the derivatives
        for(t in seq_along(hazards)) {
                derivatives[[t]] <- vector(mode = "list", length = length(lna_comp_codes))

                for(s in seq_along(derivatives[[t]])) {
                        derivatives[[t]][[s]] <- D(rate_syms[[t]], comp_syms[[s]])
                        derivatives[[t]][[s]] <- paste(deparse(derivatives[[t]][[s]]), collapse = "")
                }

                if(!is.na(time_ind)) {
                        time_derivs[[t]] <- D(rate_syms[[t]], time_sym)
                        time_derivs[[t]] <- paste(deparse(time_derivs[[t]]), collapse = "")
                } else {
                        time_derivs[[t]] <- "0.0"
                }
        }
        derivatives <- unlist(derivatives)
        time_derivs <- unlist(time_derivs)

        # # attempt to simplify the hazards and derivatives if possible
        # for(s in seq_along(hazards)) {
        #         hazards[s]     <- as.character(Ryacas::yacas(Ryacas::Simplify(hazards[s])))
        #         time_derivs[s] <- as.character(Ryacas::yacas(Ryacas::Simplify(time_derivs[s])))
        # }
        #
        # for(s in seq_along(derivatives)) {
        #         derivatives[s] <- as.character(Ryacas::yacas(Ryacas::Simplify(derivatives[s])))
        # }

        # replace the hash codes with the names of the vector elements
        for(s in seq_along(hazards)) {
                for(j in seq_len(nrow(lookup_table))) {
                        hazards[s]      <- gsub(pattern = lookup_table[j,"code"], replacement = lookup_table[j,"varname"], x = hazards[s])
                        hazards[s]      <- gsub(" ", "", hazards[s])
                        time_derivs[s]  <- gsub(pattern = lookup_table[j,"code"], replacement = lookup_table[j,"varname"], x = time_derivs[s])
                }
        }

        for(s in seq_along(derivatives)) {
                for(j in seq_len(nrow(lookup_table))) {
                        derivatives[s] <- gsub(pattern = lookup_table[j,"code"], replacement = lookup_table[j,"varname"], x = derivatives[s])
                        derivatives[s] <- gsub(" ", "", derivatives[s])
                }
                derivatives[s] <- sub_powers(derivatives[s])
        }

        # substitute exp(Z[*]), exp(-Z[*]), and exp(-2*Z[*]) for precomputed vector elements
        for(s in seq_along(hazards)) {
                # exp(Z) expressions
                exp_Z_matches <- gregexpr("exp\\(Z\\[\\d\\]\\)", hazards[s])
                exp_Z_indices <- unlist(regmatches(hazards[s], exp_Z_matches))
                exp_Z_indices <- as.character(unlist(regmatches(exp_Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_Z_indices))))

                hazards[s] <- gsub(pattern = "exp\\(Z\\[\\d\\]\\)", "exp_Z__INDEX__", hazards[s])
                for(r in seq_along(exp_Z_indices)) {
                        hazards[s] <- sub("__INDEX__", exp_Z_indices[r], hazards[s])
                }

                # exp(-Z) expressions
                exp_negZ_matches <- gregexpr("exp\\(-Z\\[\\d\\]\\)", hazards[s])
                exp_negZ_indices <- unlist(regmatches(hazards[s], exp_negZ_matches))
                exp_negZ_indices <- as.character(unlist(regmatches(exp_negZ_indices, gregexpr("\\[[[:digit:]]+\\]", exp_negZ_indices))))

                hazards[s] <- gsub(pattern = "exp\\(-Z\\[\\d\\]\\)", "exp_negZ__INDEX__", hazards[s])
                for(r in seq_along(exp_negZ_indices)) {
                        hazards[s] <- sub("__INDEX__", exp_negZ_indices[r], hazards[s])
                }

                # exp(-2*Z) expressions
                exp_neg2Z_matches <- gregexpr("exp\\(-2\\*Z\\[\\d\\]\\)", hazards[s])
                exp_neg2Z_indices <- unlist(regmatches(hazards[s], exp_neg2Z_matches))
                exp_neg2Z_indices <- as.character(unlist(regmatches(exp_neg2Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_neg2Z_indices))))

                hazards[s] <- gsub(pattern = "exp\\(-2\\*Z\\[\\d\\]\\)", "exp_neg2Z__INDEX__", hazards[s])
                for(r in seq_along(exp_neg2Z_indices)) {
                        hazards[s] <- sub("__INDEX__", exp_neg2Z_indices[r], hazards[s])
                }
        }

        for(s in seq_along(derivatives)) {
                # exp(Z) expressions
                exp_Z_matches <- gregexpr("exp\\(Z\\[\\d\\]\\)", derivatives[s])
                exp_Z_indices <- unlist(regmatches(derivatives[s], exp_Z_matches))
                exp_Z_indices <- as.character(unlist(regmatches(exp_Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_Z_indices))))

                derivatives[s] <- gsub(pattern = "exp\\(Z\\[\\d\\]\\)", "exp_Z__INDEX__", derivatives[s])
                for(r in seq_along(exp_Z_indices)) {
                        derivatives[s] <- sub("__INDEX__", exp_Z_indices[r], derivatives[s])
                }

                # exp(-Z) expressions
                exp_negZ_matches <- gregexpr("exp\\(-Z\\[\\d\\]\\)", derivatives[s])
                exp_negZ_indices <- unlist(regmatches(derivatives[s], exp_negZ_matches))
                exp_negZ_indices <- as.character(unlist(regmatches(exp_negZ_indices, gregexpr("\\[[[:digit:]]+\\]", exp_negZ_indices))))

                derivatives[s] <- gsub(pattern = "exp\\(-Z\\[\\d\\]\\)", "exp_negZ__INDEX__", derivatives[s])
                for(r in seq_along(exp_negZ_indices)) {
                        derivatives[s] <- sub("__INDEX__", exp_negZ_indices[r], derivatives[s])
                }

                # exp(-2*Z) expressions
                exp_neg2Z_matches <- gregexpr("exp\\(-2\\*Z\\[\\d\\]\\)", derivatives[s])
                exp_neg2Z_indices <- unlist(regmatches(derivatives[s], exp_neg2Z_matches))
                exp_neg2Z_indices <- as.character(unlist(regmatches(exp_neg2Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_neg2Z_indices))))

                derivatives[s] <- gsub(pattern = "exp\\(-2\\*Z\\[\\d\\]\\)", "exp_neg2Z__INDEX__", derivatives[s])
                for(r in seq_along(exp_neg2Z_indices)) {
                        derivatives[s] <- sub("__INDEX__", exp_neg2Z_indices[r], derivatives[s])
                }
        }

        return(list(hazards = hazards, derivatives = derivatives, lna_param_codes = lna_param_codes))
}