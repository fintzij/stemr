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
      
      lna_param_codes <- c(param_codes, const_codes + length(param_codes), tcovar_codes + length(param_codes) + length(const_codes) - 1)
      
      lookup_table <- 
         data.frame(varname     = c(paste("odeintr::pars[", lna_param_codes, "]", sep = ""),
                                    paste("Z[", lna_comp_codes, "]", sep = "")),
                    search_name = c(names(param_codes),
                                    names(const_codes),
                                    names(tcovar_codes),
                                    names(lna_comp_codes)),
                    code        = NA,
                    log_code    = NA,
                    stringsAsFactors = FALSE)
      
      # get indices for which rows correspond to the compartments
      if("TIME" %in% names(tcovar_codes)){
            lookup_table[which(lookup_table[,"search_name"] == "TIME"), 1] <- "t"
      }
      comp_inds <- unname(sapply(names(lna_comp_codes), match, table = lookup_table[,"search_name"]))
      time_ind  <- match("t", lookup_table[,1])
      
      # generate the code strings
      lookup_table$code <- 
         sapply(seq_len(nrow(lookup_table)),
                function(i) paste(sample(c(letters, LETTERS), 15, replace = TRUE), collapse = ""))
      
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
      
      lookup_table$log_code[comp_inds] <- paste0("(exp(", lookup_table$code[comp_inds], ")-1)")
      
      # make the substitutions in the rate strings
      for(s in seq_along(lna_rates)) {
            for(j in seq_len(nrow(lookup_table))) {
                  lna_rates[s] <- gsub(pattern = paste0('\\<', lookup_table[j,"search_name"], '\\>'),
                                       replacement = lookup_table[j,"code"], x = lna_rates[s])
            }
            
            for(j in seq_along(comp_inds)) {
                  lna_rates[s] <- gsub(pattern = paste0('\\<', lookup_table[comp_inds[j], "code"], '\\>'),
                                       replacement = lookup_table[comp_inds[j], "log_code"], x = lna_rates[s])
            }
      }
      
      # generate hazards and derivatives for the Jacobian
      hazards     <- lna_rates
      ito_coefs   <- vector(mode = "character", length = length(hazards))
      derivatives <- vector(mode = "list", length = length(hazards))
      
      for(r in seq_along(hazards)) {
            
            ito_coefs[r] <- paste0("(exp(-",
                                   lookup_table[comp_inds[r], "code"],
                                   ") - 0.5*exp(-2*",
                                   lookup_table[comp_inds[r], "code"],
                                   "))")
            
            hazards[r] <- paste0(ito_coefs[r],"*(",hazards[r],")")
            
            hazards[r] <- gsub("\\+ \\-", "- ", hazards[r])
      }
      
      # generate symbolic expressions for the rates and other objects
      rate_syms   <- lapply(hazards, function(x) (parse(text = x)))
      comp_syms   <- lapply(lookup_table[comp_inds,"code"], Ryacas::ysym)
      time_derivs <- vector(mode = "list", length = length(hazards))
      
      # if there are time-varying rates, generate the time-derivative symbols
      if(!is.na(time_ind)) time_sym  <- Ryacas::ysym(lookup_table[time_ind,"code"])
      
      # compute the derivatives
      for(t in seq_along(hazards)) {
            derivatives[[t]] <- vector(mode = "list", length = length(lna_comp_codes))
            
            for(s in seq_along(derivatives[[t]])) {
               derivatives[[t]][[s]] <- D(rate_syms[[t]], as.character(comp_syms[[s]]))
               derivatives[[t]][[s]] <- paste(deparse(derivatives[[t]][[s]]), collapse = "")
            }
            
            if(!is.na(time_ind)) {
               time_derivs[[t]] <- D(rate_syms[[t]], as.character(time_sym))
               time_derivs[[t]] <- paste(deparse(time_derivs[[t]]), collapse = "")
            } else {
                  time_derivs[[t]] <- "0.0"
            }
      }
      
      derivatives <- unlist(derivatives)
      time_derivs <- unlist(time_derivs)
      
      # replace the hash codes with the names of the vector elements
      for(s in seq_along(lna_rates)) {
            for(j in seq_len(nrow(lookup_table))) {
                  lna_rates[s] <- 
                     gsub(pattern = paste0('\\<', lookup_table[j,"code"], '\\>'),
                          replacement = lookup_table[j,"varname"], x = lna_rates[s])
            }
         
            lna_rates[s] <- 
               paste0(deparse(sub_powers(parse(text = lna_rates[s]))[[1]]), collapse = "")
            lna_rates[s] <- gsub(" ", "", lna_rates[s])
      }
      
      for(s in seq_along(ito_coefs)) {
            for(j in seq_len(nrow(lookup_table))) {
                  ito_coefs[s] <- 
                     gsub(pattern = paste0('\\<', lookup_table[j,"code"], '\\>'),
                          replacement = lookup_table[j,"varname"], x = ito_coefs[s])
            }
         
         ito_coefs[s] <- 
            paste0(deparse(sub_powers(parse(text = ito_coefs[s]))[[1]]), collapse = "")
         ito_coefs[s] <- gsub(" ", "", ito_coefs[s])
      }
      
      for(s in seq_along(hazards)) {
            for(j in seq_len(nrow(lookup_table))) {
                  hazards[s] <-
                     gsub(pattern = paste0('\\<', lookup_table[j,"code"], '\\>'),
                          replacement = lookup_table[j,"varname"], x = hazards[s])
                  
                  time_derivs[s] <-
                     gsub(pattern = paste0('\\<', lookup_table[j,"code"], '\\>'), 
                          replacement = lookup_table[j,"varname"], x = time_derivs[s])
            }
         
         hazards[s] <- 
            paste0(deparse(sub_powers(parse(text = hazards[s]))[[1]]), collapse = "")
         hazards[s] <- gsub(" ", "", hazards[s])
         
         time_derivs[s] <- 
            paste0(deparse(sub_powers(parse(text = time_derivs[s]))[[1]]), collapse = "")
         time_derivs[s] <- gsub(" ", "", time_derivs[s])
      }
      
      for(s in seq_along(derivatives)) {
            for(j in seq_len(nrow(lookup_table))) {
                  derivatives[s] <- 
                     gsub(pattern = paste0('\\<', lookup_table[j,"code"], '\\>'), 
                          replacement = lookup_table[j,"varname"], x = derivatives[s])
            }
         derivatives[s] <- 
            paste0(deparse(sub_powers(parse(text = derivatives[s]))[[1]]), collapse = "")
         derivatives[s] <- gsub(" ", "", derivatives[s])
      }
      
      # substitute exp(Z[*]), expm1(Z[*]), exp(-Z[*]), and exp(-2*Z[*]) for precomputed vector elements
      for(s in seq_along(ito_coefs)) {
            # exp(-Z) expressions
            exp_neg_Z_matches <- gregexpr("exp\\(-Z\\[[[:digit:]]+\\]\\)", ito_coefs[s])
            exp_neg_Z_indices <- unlist(regmatches(ito_coefs[s], exp_neg_Z_matches))
            exp_neg_Z_indices <- as.character(unlist(regmatches(exp_neg_Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_neg_Z_indices))))
            
            ito_coefs[s] <- gsub(pattern = "exp\\(-Z\\[[[:digit:]]+\\]\\)", "odeintr::exp_neg_Z__INDEX__", ito_coefs[s])
            for(r in seq_along(exp_neg_Z_indices)) {
                  ito_coefs[s] <- sub("__INDEX__", exp_neg_Z_indices[r], ito_coefs[s])
            }
            
            # exp(-2*Z) expressions
            exp_neg_2Z_matches <- gregexpr("exp\\(-2\\*Z\\[[[:digit:]]+\\]\\)", ito_coefs[s])
            exp_neg_2Z_indices <- unlist(regmatches(ito_coefs[s], exp_neg_2Z_matches))
            exp_neg_2Z_indices <- as.character(unlist(regmatches(exp_neg_2Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_neg_2Z_indices))))
            
            ito_coefs[s] <- gsub(pattern = "exp\\(-2\\*Z\\[[[:digit:]]+\\]\\)", "odeintr::exp_neg_2Z__INDEX__", ito_coefs[s])
            for(r in seq_along(exp_neg_2Z_indices)) {
                  ito_coefs[s] <- sub("__INDEX__", exp_neg_2Z_indices[r], ito_coefs[s])
            }
      }
      
      for(s in seq_along(lna_rates)) {
            
            # expm1(Z) expressions
            expm1_Z_matches <- gregexpr("\\(exp\\(Z\\[[[:digit:]]+\\]\\)-1\\)", lna_rates[s])
            expm1_Z_indices <- unlist(regmatches(lna_rates[s], expm1_Z_matches))
            expm1_Z_indices <- as.character(unlist(regmatches(expm1_Z_indices, gregexpr("\\[[[:digit:]]+\\]", expm1_Z_indices))))
            
            lna_rates[s] <- gsub(pattern = "\\(exp\\(Z\\[[[:digit:]]+\\]\\)-1\\)", "odeintr::expm1_Z__INDEX__", lna_rates[s])
            for(r in seq_along(expm1_Z_indices)) {
                  lna_rates[s] <- sub("__INDEX__", expm1_Z_indices[r], lna_rates[s])
            }
            
            # exp(Z) expressions
            exp_Z_matches <- gregexpr("exp\\(Z\\[[[:digit:]]+\\]\\)", lna_rates[s])
            exp_Z_indices <- unlist(regmatches(lna_rates[s], exp_Z_matches))
            exp_Z_indices <- as.character(unlist(regmatches(exp_Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_Z_indices))))
            
            lna_rates[s] <- gsub(pattern = "exp\\(Z\\[[[:digit:]]+\\]\\)", "odeintr::exp_Z__INDEX__", lna_rates[s])
            for(r in seq_along(exp_Z_indices)) {
                  lna_rates[s] <- sub("__INDEX__", exp_Z_indices[r], lna_rates[s])
            }
            
            # exp(-Z) expressions
            exp_neg_Z_matches <- gregexpr("exp\\(-Z\\[[[:digit:]]+\\]\\)", lna_rates[s])
            exp_neg_Z_indices <- unlist(regmatches(lna_rates[s], exp_neg_Z_matches))
            exp_neg_Z_indices <- as.character(unlist(regmatches(exp_neg_Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_neg_Z_indices))))
            
            lna_rates[s] <- gsub(pattern = "exp\\(-Z\\[[[:digit:]]+\\]\\)", "odeintr::exp_neg_Z__INDEX__", lna_rates[s])
            for(r in seq_along(exp_neg_Z_indices)) {
                  lna_rates[s] <- sub("__INDEX__", exp_neg_Z_indices[r], lna_rates[s])
            }
            
            # exp(-2*Z) expressions
            exp_neg_2Z_matches <- gregexpr("exp\\(-2\\*Z\\[[[:digit:]]+\\]\\)", lna_rates[s])
            exp_neg_2Z_indices <- unlist(regmatches(lna_rates[s], exp_neg_2Z_matches))
            exp_neg_2Z_indices <- as.character(unlist(regmatches(exp_neg_2Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_neg_2Z_indices))))
            
            lna_rates[s] <- gsub(pattern = "exp\\(-2\\*Z\\[[[:digit:]]+\\]\\)", "odeintr::exp_neg_2Z__INDEX__", lna_rates[s])
            for(r in seq_along(exp_neg_2Z_indices)) {
                  lna_rates[s] <- sub("__INDEX__", exp_neg_2Z_indices[r], lna_rates[s])
            }
      }
      
      for(s in seq_along(hazards)) {
            
            # expm1(Z) expm1ressions
            expm1_Z_matches <- gregexpr("\\(exp\\(Z\\[[[:digit:]]+\\]\\)-1\\)", hazards[s])
            expm1_Z_indices <- unlist(regmatches(hazards[s], expm1_Z_matches))
            expm1_Z_indices <- as.character(unlist(regmatches(expm1_Z_indices, gregexpr("\\[[[:digit:]]+\\]", expm1_Z_indices))))
            
            hazards[s] <- gsub(pattern = "\\(exp\\(Z\\[[[:digit:]]+\\]\\)-1\\)", "odeintr::expm1_Z__INDEX__", hazards[s])
            for(r in seq_along(expm1_Z_indices)) {
                  hazards[s] <- sub("__INDEX__", expm1_Z_indices[r], hazards[s])
            }
            
            # exp(Z) expressions
            exp_Z_matches <- gregexpr("exp\\(Z\\[[[:digit:]]+\\]\\)", hazards[s])
            exp_Z_indices <- unlist(regmatches(hazards[s], exp_Z_matches))
            exp_Z_indices <- as.character(unlist(regmatches(exp_Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_Z_indices))))
            
            hazards[s] <- gsub(pattern = "exp\\(Z\\[[[:digit:]]+\\]\\)", "odeintr::exp_Z__INDEX__", hazards[s])
            for(r in seq_along(exp_Z_indices)) {
                  hazards[s] <- sub("__INDEX__", exp_Z_indices[r], hazards[s])
            }
            
            # exp(-Z) expressions
            exp_neg_Z_matches <- gregexpr("exp\\(-Z\\[[[:digit:]]+\\]\\)", hazards[s])
            exp_neg_Z_indices <- unlist(regmatches(hazards[s], exp_neg_Z_matches))
            exp_neg_Z_indices <- as.character(unlist(regmatches(exp_neg_Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_neg_Z_indices))))
            
            hazards[s] <- gsub(pattern = "exp\\(-Z\\[[[:digit:]]+\\]\\)", "odeintr::exp_neg_Z__INDEX__", hazards[s])
            for(r in seq_along(exp_neg_Z_indices)) {
                  hazards[s] <- sub("__INDEX__", exp_neg_Z_indices[r], hazards[s])
            }
            
            # exp(-2*Z) expressions
            exp_neg_2Z_matches <- gregexpr("exp\\(-2\\*Z\\[[[:digit:]]+\\]\\)", hazards[s])
            exp_neg_2Z_indices <- unlist(regmatches(hazards[s], exp_neg_2Z_matches))
            exp_neg_2Z_indices <- as.character(unlist(regmatches(exp_neg_2Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_neg_2Z_indices))))
            
            hazards[s] <- gsub(pattern = "exp\\(-2\\*Z\\[[[:digit:]]+\\]\\)", "odeintr::exp_neg_2Z__INDEX__", hazards[s])
            for(r in seq_along(exp_neg_2Z_indices)) {
                  hazards[s] <- sub("__INDEX__", exp_neg_2Z_indices[r], hazards[s])
            }
      }
      
      for(s in seq_along(derivatives)) {
            
            # exp(Z) expressions
            expm1Z_matches <- gregexpr("\\(exp\\(Z\\[[[:digit:]]+\\]\\)-1\\)", derivatives[s])
            expm1Z_indices <- unlist(regmatches(derivatives[s], expm1Z_matches))
            expm1Z_indices <- as.character(unlist(regmatches(expm1Z_indices, gregexpr("\\[[[:digit:]]+\\]", expm1Z_indices))))
            
            derivatives[s] <- gsub(pattern = "\\(exp\\(Z\\[[[:digit:]]+\\]\\)-1\\)", "odeintr::expm1_Z__INDEX__", derivatives[s])
            for(r in seq_along(expm1Z_indices)) {
                  derivatives[s] <- sub("__INDEX__", expm1Z_indices[r], derivatives[s])
            }
            
            # exp(Z) expressions
            exp_Z_matches <- gregexpr("exp\\(Z\\[[[:digit:]]+\\]\\)", derivatives[s])
            exp_Z_indices <- unlist(regmatches(derivatives[s], exp_Z_matches))
            exp_Z_indices <- as.character(unlist(regmatches(exp_Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_Z_indices))))
            
            derivatives[s] <- gsub(pattern = "exp\\(Z\\[[[:digit:]]+\\]\\)", "odeintr::exp_Z__INDEX__", derivatives[s])
            for(r in seq_along(exp_Z_indices)) {
                  derivatives[s] <- sub("__INDEX__", exp_Z_indices[r], derivatives[s])
            }
            
            # exp(-Z) expressions
            exp_neg_Z_matches <- gregexpr("exp\\(-Z\\[[[:digit:]]+\\]\\)", derivatives[s])
            exp_neg_Z_indices <- unlist(regmatches(derivatives[s], exp_neg_Z_matches))
            exp_neg_Z_indices <- as.character(unlist(regmatches(exp_neg_Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_neg_Z_indices))))
            
            derivatives[s] <- gsub(pattern = "exp\\(-Z\\[[[:digit:]]+\\]\\)", "odeintr::exp_neg_Z__INDEX__", derivatives[s])
            for(r in seq_along(exp_neg_Z_indices)) {
                  derivatives[s] <- sub("__INDEX__", exp_neg_Z_indices[r], derivatives[s])
            }
            
            # exp(-2*Z) expressions
            exp_neg_2Z_matches <- gregexpr("exp\\(-2\\*Z\\[[[:digit:]]+\\]\\)", derivatives[s])
            exp_neg_2Z_indices <- unlist(regmatches(derivatives[s], exp_neg_2Z_matches))
            exp_neg_2Z_indices <- as.character(unlist(regmatches(exp_neg_2Z_indices, gregexpr("\\[[[:digit:]]+\\]", exp_neg_2Z_indices))))
            
            derivatives[s] <- gsub(pattern = "exp\\(-2\\*Z\\[[[:digit:]]+\\]\\)", "odeintr::exp_neg_2Z__INDEX__", derivatives[s])
            for(r in seq_along(exp_neg_2Z_indices)) {
                  derivatives[s] <- sub("__INDEX__", exp_neg_2Z_indices[r], derivatives[s])
            }
      }
      
      return(list(lna_rates       = lna_rates,
                  ito_coefs       = ito_coefs,
                  hazards         = hazards,
                  derivatives     = derivatives,
                  lna_param_codes = lna_param_codes))
}
