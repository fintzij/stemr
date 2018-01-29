#' Generate objects for setting up a random walk of order 3
#'
#' @param ntimes number of times (nodes)
#'
#' @return list containing the difference matrix, the structure matrix, the
#'   singular value decomposition of the structure matrix, and the reference
#'   standard deviation of the random walk (i.e., the geometric mean standard
#'   deviation of the random walk with precision 1).
#' @export
generate_rw3 <- function(ntimes) {
      
      # Difference matrix
      D <- diag(-1, nrow = ntimes-3, ncol = ntimes)
      D[matrix(c(1:(ntimes-3), 2:(ntimes-2)), ncol = 2)] <- 3
      D[matrix(c(1:(ntimes-3), 3:(ntimes-1)), ncol = 2)] <- -3
      D[matrix(c(1:(ntimes-3), 4:ntimes), ncol = 2)] <- 1
      
      # Structure matrix
      R <- t(D)%*%D
      
      # SVD of the structure matrix
      R_svd <- svd(R)
      R_svd$d[(ntimes-1):ntimes] <- 0.0
      
      # kernel
      kern <- matrix(c(rep(1,ntimes), 1:ntimes, c(1:ntimes)^2), ncol = 3)
      kern_outer <- kern %*% t(kern)
      kern_svd <- svd(kern_outer)
      
      # reference standard deviation
      sigma_ref <- exp(sum(0.5 * log(diag(MASS::ginv(R))))/ntimes)
      
      return(list(D = D, R = R, R_svd = R_svd, kern = kern, kern_outer = kern_outer, kern_svd = kern_svd, sigma_ref = sigma_ref))
}
