#' Construct an indicator matrix for which measurement process variables are
#' measured at which observation times.
#'
#' @param obsmat observation matrix with NA indicating no observation at a
#'   particular time.
#'
#' @return logical matrix indicating measurment status
#' @export
build_measproc_indmat <- function(obsmat) {

        indmat <- !is.na(obsmat[,-1,drop=FALSE])
        colnames(indmat) <- colnames(obsmat)[-1]

        return(indmat)
}