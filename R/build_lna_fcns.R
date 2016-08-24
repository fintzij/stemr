#' Reformat the rates for use in the LNA and compute the derivatives of the rate
#' functions.
#'
#' @param rate_fcns list of rates, as returned by stem_dynamics
#' @param param_codes vector of parameter codes
#' @param const_codes vector of time-varying covariate codes
#' @param tcovar_codes vector of constant codes
#' @param comp_codes vector of compartment codes
#' @param flow_matrix matrix indicating compartment flow, whose transpose is the
#'   stoichiometry matrix
#' @param lna_scale either "log" or "linear"
#'
#' @return strings for computing the LNA hazards and the LNA jacobian
#' @export
build_lna_fcns <- function(rate_fcns, param_codes, const_codes, tcovar_codes, compartment_codes, flow_matrix, lna_scale) {

        lna_param_codes <- c(param_codes, const_codes + length(param_codes), tcovar_codes + length(param_codes) + length(const_codes)-1)
        stoich_mat      <- t(flow_matrix)
        log_scale       <- lna_scale == "log"

        lookup_table <- data.frame(varname     = c(paste("parms[", lna_param_codes, "]", sep = ""),
                                                   paste("drift[", compartment_codes, "]", sep = "")),
                                   search_name = c(paste("parameters\\[", param_codes, "\\]", sep = ""),
                                                   paste("constants\\[", const_codes, "\\]", sep = ""),
                                                   paste("tcovar\\[", tcovar_codes, "\\]", sep = ""),
                                                   paste("state\\[", compartment_codes, "\\]", sep = "")),
                                   code        = NA,
                                   log_code    = NA,
                                   stringsAsFactors = FALSE)

        # get indices for which rows correspond to the compartments
        if("TIME" %in% names(tcovar_codes)){
                lookup_table[which(lookup_table[,"search_name"] == paste("tcovar\\[",which(names(tcovar_codes) == "TIME"),"\\]", sep = "")), 1] <- "t"
        }
        comp_inds <- unname(sapply(paste("state\\[",0:(length(compartment_codes)-1), "\\]", sep = ""), match, table = lookup_table[,"search_name"]))
        time_ind  <- match("t", lookup_table[,1])

        # generate the code strings
        lookup_table$code <- lookup_table$log_code <- replicate(nrow(lookup_table),
                                                                paste(sample(c(letters, LETTERS), 15, replace = TRUE), collapse = ""),
                                                                simplify = T)
        lookup_table$log_code[comp_inds] <- paste0("exp(", lookup_table$log_code[comp_inds], ")")

        # get the unparsed rate strings
        rates <- sapply(rate_fcns, "[[", "lumped")

        # make the substitutions from the lookup table
        for(s in seq_along(rates)) {
                for(j in seq_len(nrow(lookup_table))) {
                        if(log_scale) {
                                rates[s] <- gsub(pattern = lookup_table[j,"search_name"],
                                                        replacement = lookup_table[j,"log_code"], x = rates[s])
                        } else if(!log_scale) {
                                rates[s] <- gsub(pattern = lookup_table[j,"search_name"],
                                                        replacement = lookup_table[j,"code"], x = rates[s])
                        }
                }
        }

        # generate hazards and derivatives for the Jacobian
        hazards     <- vector(mode = "character", length = length(compartment_codes))
        derivatives <- vector(mode = "list", length = length(hazards))

        if(lna_scale == "linear") {

                for(r in seq_along(hazards)) {
                        rate_inds  <- stoich_mat[r,] != 0                             # identify the rates in the linear combination
                        signs      <- ifelse(stoich_mat[r,rate_inds] == 1, "", "-")   # get the signs of the rates
                        hazards[r] <- paste(paste0(signs, rates[rate_inds]), collapse = " + ", sep = "")
                        hazards[r] <- gsub("\\+ \\-", "- ", hazards[r])
                }

        } else if(lna_scale == "log") {

                for(r in seq_along(hazards)) {
                        rate_inds  <- stoich_mat[r,] != 0                             # identify the rates in the linear combination
                        signs      <- ifelse(stoich_mat[r,rate_inds] == 1, "", "-")   # get the signs of the rates
                        hazards[r] <- paste(paste0("exp(-", lookup_table[comp_inds[r], "code"],")*("),
                                            paste0(signs, rates[rate_inds], collapse = " + ", sep = ""), ")", sep = "")
                        hazards[r] <- gsub("\\+ \\-", "- ", hazards[r])

                        # get the coefficients for the extra terms from Ito's formula
                        ito_coefs  <- stoich_mat[r,]!=0
                        hazards[r] <- paste(hazards[r], paste0(" - 0.5*exp(-2*", lookup_table[comp_inds[r], "code"],")*("),
                                            paste0(rates[ito_coefs], collapse = " + ", sep = ""), ")", sep = "")
                        hazards[r] <- as.character(Ryacas::yacas(Ryacas::Simplify(hazards[r])))

                }
        }

        # generate symbolic expressions for the rates and other objects
        rate_syms   <- lapply(hazards, function(x) (parse(text = x)))
        comp_syms   <- lapply(lookup_table[comp_inds,"code"], Ryacas::Sym)
        time_derivs <- vector(mode = "list", length = length(hazards))

        # if there are time-varying rates, generate the time-derivative symbols
        if(!is.na(time_ind)) time_sym  <- Ryacas::Sym(lookup_table[time_ind,"code"])

        # compute the derivatives
        for(t in seq_along(hazards)) {
                derivatives[[t]] <- vector(mode = "list", length = length(compartment_codes))

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
        for(s in seq_along(hazards)) {
                hazards[s]     <- as.character(Ryacas::yacas(Ryacas::Simplify(hazards[s])))
                time_derivs[s] <- as.character(Ryacas::yacas(Ryacas::Simplify(time_derivs[s])))
        }

        for(s in seq_along(derivatives)) {
                derivatives[s] <- as.character(Ryacas::yacas(Ryacas::Simplify(derivatives[s])))
        }

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
        }

        for(s in seq_along(rates)) {
                for(j in seq_len(nrow(lookup_table))) {
                        rates[s]      <- gsub(pattern = lookup_table[j,"code"], replacement = lookup_table[j,"varname"], x = rates[s])
                        rates[s]      <- gsub(" ", "", rates[s])
                }
        }

        return(list(rates = rates, hazards = hazards, derivatives = derivatives, time_derivs = time_derivs, lna_param_codes = lna_param_codes))
}