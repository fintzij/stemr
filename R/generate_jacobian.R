generate_jacobian <- function(rates, parameters, tcovar, constants, compartments) {

        # Underscores will need to be removed and re-inserted at the end.
        # Thus, we create a lookup table for the strings with and without underscores
        lookup_table      <- data.frame(c(parameters, tcovar, constants, compartments), code = NA)
        lookup_table      <- lookup_table[order(sapply(as.character(lookup_table[,1]), nchar), decreasing = T),]
        lookup_table$code <- replicate(nrow(lookup_table), paste(sample(c(letters, LETTERS), 15, replace = TRUE), collapse = ""), simplify = T)

        # get indices for which rows correspond to the compartments
        comp_inds <- sapply(compartments, match, table = lookup_table[,1])

        # get the unparsed rate strings
        rate_strings <- sapply(rates, "[[", "unparsed")

        # make the substitutions from the lookup table
        for(s in seq_along(rate_strings)) {
                for(j in seq_len(nrow(lookup_table))) {
                        rate_strings[s] <- gsub(pattern = lookup_table[j,1], replacement = lookup_table[j,2], x = rate_strings[s])
                }
        }

        # generate symbolic expressions for the rates and other objects
        rate_syms   <- lapply(rate_strings, Ryacas::Sym)
        object_syms <- lapply(lookup_table[comp_inds,2], Ryacas::Sym)

        # create list to store derivatives
        derivatives <- vector(mode = "list", length = length(rate_strings))

        # compute the derivatives
        for(t in seq_along(rate_strings)) {
                derivatives[[t]] <- lapply(object_syms, FUN = function(x) Ryacas::deriv.Sym(rate_syms[[t]], x))
        }

        derivatives   <- unlist(derivatives, recursive = FALSE)

        deriv_strings <- sapply(derivatives)
}