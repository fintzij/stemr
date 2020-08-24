#' Parse a string and substitute powers of the form a^b with pow(a,b).
#'
#' The function for accomplishing this was proposed by Roland on the following
#' \href{https://stackoverflow.com/questions/40606723/substitute-the-power-symbol-with-cs-pow-syntax-in-mathematical-expression/40612497}{stack
#' overflow thread.}
#'
#' @param string
#'
#' @return string suitable for C++
#' @export
sub_powers <- function(string) {    
        
        #check if you are at the end of the tree's branch
        if (is.name(string) || is.atomic(string)) { 
                
                #replace ^
                if (string == quote(`^`)) return(quote(pow))
                
                return(string)
        }
        
        if (string[[1]] == quote(sqrt)) {
                
                #replace sqrt
                string[[1]] <- quote(pow)
                
                #add the second argument
                string[[3]] <- quote(0.5)
        }
        
        #follow the tree with recursion
        for (i in seq_along(string)) string[[i]] <- sub_powers(string[[i]])
        
        return(string)    
}
