#' Parse a string and substitute powers of the form a^b with pow(a,b).
#'
#' @param string
#'
#' @return string suitable for C++
#' @export
sub_powers <- function(string) {

        # remove any white spaces
        string <- gsub(" ", "", string)

        while(grepl(pattern = "\\^", x = string)) {

                # get the number of powers to replace
                matches <- gregexpr("\\^", string)

                # get the substring for making replacements
                psym <- matches[[1]][[1]] # location of the first ^

                # get the power
                m_start <- psym + 1
                n_open <- 0; len <- 0
                while(n_open >= 0) {
                        if(substr(string, m_start + len, m_start + len) == "(") {
                                n_open <- n_open + 1
                        } else if(substr(string, m_start + len, m_start + len) == ")") {
                                n_open <- n_open - 1
                        }
                      
                        if(n_open != 0) {
                                len <- len + 1
                        } else if(n_open == 0) {
                                m_end <- m_start + len
                        }
                }

                m_start_pow <- m_start+1; m_end_pow <- m_end-1
                m_fcn_pow <- unlist(strsplit(gsub(" ", "", substr(string, m_start_pow, m_end_pow)), ","))

                # get the base
                m_start <- psym - 1
                n_open <- 0; len <- 0
                while(n_open >= 0) {
                        if(substr(string, m_start - len, m_start - len) == ")") {
                                n_open <- n_open + 1
                        } else if(substr(string, m_start - len, m_start - len) == "(") {
                                n_open <- n_open - 1
                        }
                        if(n_open != 0) {
                                len <- len + 1
                        } else if(n_open == 0) {
                                m_end <- m_start - len
                        }
                }

                m_start_base <- m_end+1; m_end_base <- m_start-1
                m_fcn_base <- unlist(strsplit(gsub(" ", "", substr(string, m_start_base, m_end_base)), ","))

                # construct and make the replacement
                replacement <- paste0("pow(",m_fcn_base, ",", m_fcn_pow,")")
                string      <- sub(pattern = substr(string, m_start_base-1, m_end_pow+1),
                                    replacement = replacement, x = string, fixed = TRUE)
        }

        return(string)
}