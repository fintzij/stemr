#' Parse the ODE rates so they can be compiled.
#'
#' @param ode_rates character vector of ODE rates
#' @param param_codes named numeric vector of parameter codes
#' @param const_codes named numeric vector of constant codes
#' @param tcovar_codes named numeric vector of time-varying covariate codes
#' @param ode_comp_codes named numeric vector of ODE compartment codes
#'
#' @return string snippets for the ODE that can be compiled
#' @export
parse_ode_rates <- function(ode_rates, param_codes, const_codes, tcovar_codes, ode_comp_codes) {

        ode_param_codes <- c(param_codes, const_codes + length(param_codes), tcovar_codes + length(param_codes) + length(const_codes) - 1)

        lookup_table <- data.frame(varname     = c(paste("odeintr::pars[", ode_param_codes, "]", sep = ""),
                                                   paste("x[", ode_comp_codes, "]", sep = "")),
                                   search_name = c(names(param_codes),
                                                   names(const_codes),
                                                   names(tcovar_codes),
                                                   names(ode_comp_codes)),
                                   code        = NA,
                                   stringsAsFactors = FALSE)

        # get indices for which rows correspond to the compartments
        if("TIME" %in% names(tcovar_codes)){
                lookup_table[which(lookup_table[,"search_name"] == "TIME"), 1] <- "t"
        }
        comp_inds <- unname(sapply(names(ode_comp_codes), match, table = lookup_table[,"search_name"]))
        time_ind  <- match("t", lookup_table[,1])

        # generate the code strings
        lookup_table$code <- 
                sapply(seq_len(nrow(lookup_table)),
                       function(i) paste(sample(c(letters, LETTERS), 15, replace = TRUE), 
                                         collapse = ""))

        # make sure there are no partial matches between columns in the lookup table
        while(any(sapply(lookup_table[,"search_name"], grepl, x = lookup_table[,"code"]))) {
                which_match <- 
                        which(apply(sapply(lookup_table[,"search_name"], 
                                           grepl, x = lookup_table[,"code"]), 1, any))
                
                for(m in seq_along(which_match)) {
                        lookup_table[which_match[m],"code"] <- 
                                paste(sample(c(letters,LETTERS), 15, replace = TRUE), collapse = "")
                }
        }
        
        while(any(duplicated(lookup_table$code))) {
                lookup_table$code[duplicated(lookup_table$code)] = 
                        sapply(seq_len(sum(duplicated(lookup_table$code))),
                               function(i) paste(sample(c(letters, LETTERS), 15, replace = TRUE), collapse = ""))
        }

        # make the substitutions in the rate strings
        for(s in seq_along(ode_rates)) {
                for(j in seq_len(nrow(lookup_table))) {
                        ode_rates[s] <- gsub(pattern = paste0('\\<', lookup_table[j,"search_name"], '\\>'),
                                             replacement = lookup_table[j,"code"], x = ode_rates[s])
                }
        }

        # generate hazards and derivatives for the Jacobian
        hazards     <- ode_rates

        # replace the hash codes with the names of the vector elements
        for(s in seq_along(hazards)) {
                for(j in seq_len(nrow(lookup_table))) {
                        hazards[s] <- 
                                gsub(pattern = paste0('\\<',lookup_table[j,"code"],'\\>'),
                                     replacement = lookup_table[j,"varname"], x = hazards[s])
                }
                
                hazards[s] <- 
                        paste0(deparse(sub_powers(parse(text = hazards[s]))[[1]]), collapse = "")
                hazards[s] <- gsub(" ", "", hazards[s])
        }

        for(s in seq_along(ode_rates)) {
                for(j in seq_len(nrow(lookup_table))) {
                        ode_rates[s] <- 
                                gsub(pattern = paste0('\\<',lookup_table[j,"code"],'\\>'),
                                     replacement = lookup_table[j,"varname"], x = ode_rates[s])
                }
                
                ode_rates[s] <- 
                        paste0(deparse(sub_powers(parse(text = ode_rates[s]))[[1]]), collapse = "")
                ode_rates[s] <- gsub(" ", "", ode_rates[s])
        }

        return(list(hazards = hazards, ode_param_codes = ode_param_codes))
}
