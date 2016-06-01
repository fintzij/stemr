#' Makes comp_fcn replacements.
#'
#' @param rate_string rate string in which to make replacements
#' @param comp character string to be replaced
#' @param subs character vector of compartments to be substituted
#'
#' @return rate string
#' @export
sub_comp_fcns <- function(rate_string, comp, subs) {

        # get the number of comp_strings to replace
        matches <- gregexpr("comp_fcn", rate_string)

        for(m in seq_along(matches)) {

                # get the substring for making replacements
                m_start <- matches[[m]][[1]] + 9
                n_open <- 1; len <- 0
                while(n_open > 0) {
                        if(substr(rate_string, m_start + len, m_start + len) == "(") {
                                n_open <- n_open + 1
                        } else if(substr(rate_string, m_start + len, m_start + len) == ")") {
                                n_open <- n_open - 1
                        }
                        if(n_open != 0) {
                                len <- len + 1
                        } else if(n_open == 0) {
                                m_end <- m_start + len - 1
                        }
                }

                m_fcn <- unlist(strsplit(gsub(" ", "", substr(rate_string, m_start, m_end)), ","))

                # make substitutions
                rep_fcn <- rep(m_fcn[1], length(subs))
                for(s in seq_along(rep_fcn)) {
                        rep_fcn[s] <- gsub(comp, subs[s], rep_fcn[s])
                }

                # combine
                if(m_fcn[2] == "sum") {
                        sub_fcn <- paste0("(", paste(rep_fcn, collapse = " + "), ")")
                } else if(m_fcn[2] == "prod") {
                        sub_fcn <- paste0("(", paste(rep_fcn, collapse = " * "), ")")
                } else {
                        # check that aggregation is either sum or prod
                        comp_fcn(sub_fcn, m_fcn[2])
                }

                # replace
                rate_string <- gsub(pattern = substr(rate_string, matches[[m]][[1]], m_end + 1), replacement = sub_fcn, x =  rate_string, fixed = TRUE)
        }

        return(rate_string)
}