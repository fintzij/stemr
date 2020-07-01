#' Wrapper for evaluating the compartment counts at census times.
#'
#' @param path matrix containing the path to be censused.
#' @param census_times vector of census times
#' @param census_columns indices of columns to be censused. If not specified,
#'   will default to all columns except 'time', 'event', and 'ID'.
#' @param census_matrix existing matrix of compartment counts at census times to
#'   be modified in place. If not specified, a new matrix will be returned.
#'
#' @return matrix of compartment counts at census times
#' @export
census_path <- function(path, census_times, census_columns = NULL, census_matrix = NULL) {

        if(is.unsorted(path[,1])) {
                if(is.null(census_matrix)) {
                        warning("The supplied path matrix had to be sorted. It would be faster to provide an already sorted path matrix.")
                        path <- path[order(path[,1]), , drop = FALSE]
                } else {
                        stop("The path matrix must be pre-sorted if a census matrix is to be modified in place.")
                }
        }

        # get the names for columns to be censused.
        path_colnames <- colnames(path)
        if(is.null(census_columns)) {
                census_columns  <- which(!path_colnames %in% c("time", "event", "ID"))
        }
        census_colnames <- path_colnames[census_columns]

        if(is.null(census_matrix)) {
                census_matrix <- build_census_path(path, census_times, census_columns - 1)
                colnames(census_matrix) <- c("time", census_colnames)
                return(census_matrix)

        } else {
                retrieve_census_path(census_matrix, path, census_times, census_columns - 1)
        }
}