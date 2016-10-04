#' Parse the LNA rates so they can be compiled.
#'
#' @param lna_rates character vector of LNA rates
#' @param param_codes named numeric vector of parameter codes
#' @param const_codes named numeric vector of constant codes
#' @param tcovar_codes named numeric vector of time-varying covariate codes
#' @param lna_comp_codes named numeric vector of LNA compartment codes
#'
#' @return string snippets for the LNA that can be compiled
#' @export
parse_lna_rates <- function(lna_rates, param_codes, const_codes, tcovar_codes, lna_comp_codes) {

        lna_param_codes <- c(param_codes, const_codes + length(param_codes), tcovar_codes + length(param_codes) + length(const_codes)-1)

        lookup_table <- data.frame(varname     = c(paste("pars[", lna_param_codes, "]", sep = ""),
                                                   paste("x[", lna_comp_codes, "]", sep = "")),
                                   search_name = c(names(param_codes),
                                                   names(const_codes),
                                                   names(tcovar_codes),
                                                   names(lna_comp_codes)),
                                   code        = NA,
                                   stringsAsFactors = FALSE)

        # get indices for which rows correspond to the compartments
        if("TIME" %in% names(tcovar_codes)){
                lookup_table[which(lookup_table[,"search_name"] == "TIME"), 1] <- "t"
        }
        comp_inds <- unname(sapply(names(lna_comp_codes), match, table = lookup_table[,"search_name"]))
        time_ind  <- match("t", lookup_table[,1])

        # generate the code strings
        lookup_table$code <- replicate(nrow(lookup_table),
                                       paste(sample(c(letters, LETTERS), 15, replace = TRUE), collapse = ""),
                                       simplify = T)


        # make the substitutions in the rate strings
        for(s in seq_along(lna_rates)) {
                for(j in seq_len(nrow(lookup_table))) {
                        lna_rates[s] <- gsub(pattern = lookup_table[j,"search_name"],
                                             replacement = lookup_table[j,"code"], x = lna_rates[s])
                }
        }

        # generate hazards and derivatives for the Jacobian
        hazards     <- lna_rates
        derivatives <- vector(mode = "list", length = length(hazards))

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

        # attempt to simplify the hazards and derivatives if possible
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

        for(s in seq_along(lna_rates)) {
                for(j in seq_len(nrow(lookup_table))) {
                        lna_rates[s]      <- gsub(pattern = lookup_table[j,"code"],
                                                  replacement = lookup_table[j,"varname"], x = lna_rates[s])
                        lna_rates[s]      <- gsub(" ", "", lna_rates[s])
                }
                lna_rates[s] <- sub_powers(lna_rates[s])
        }

        return(list(rates = lna_rates, derivatives = derivatives, lna_param_codes = lna_param_codes))
}