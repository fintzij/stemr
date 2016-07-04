#' Census each path in a collection of paths to obtain the compartment counts at
#' census times.
#'
#' @param paths either a list or array of paths
#' @inheritParams census_path
#' @param as_array return the matrices of compartment counts as an array?
#'
#' @return a list (default) or array of matrices of compartment counts at census
#'   times.
#' @export
census_path_collection <- function(paths, census_times, census_columns = NULL, as_array = FALSE) {

        if(is.list(paths)) {
                census_paths <- lapply(X = paths, FUN = census_path, census_times = census_times, census_columns = census_columns)
        } else {
                census_paths <- vector(mode = "list", length = dim(paths)[3])
                for(k in seq_len(dim(paths)[3])) {
                        census_paths[[k]] <- census_path(paths[,,k], census_times = census_times, census_columns = census_columns)
                }
        }

        if(as_array) census_paths <- array(unlist(census_paths), dim = c(length(census_times), dim(census_paths[[1]])[2], dim(paths)[3]))

        return(census_paths)

}