#' Generate objects for setting up a random walk of order 1
#'
#' @param ntimes number of times (nodes)
#'
#' @return list containing the difference matrix, the structure matrix, the
#'   singular value decomposition of the structure matrix, and the reference
#'   standard deviation of the random walk (i.e., the geometric mean standard
#'   deviation of the random walk with precision 1).
#' @export
generate_rw1 <- function(ntimes) {
      
      # Difference matrix
      D <- diag(-1, nrow = ntimes-1, ncol = ntimes)
      D[matrix(c(1:(ntimes-1), 2:(ntimes)), ncol = 2)] <- 1
      
      # Structure matrix
      R <- t(D)%*%D
      
      # SVD of the structure matrix
      R_svd <- svd(R)
      R_svd$d[ntimes] <- 0.0
      
      # kernel
      kern <- matrix(rep(1,ntimes), ncol = 1)
      kern_outer <- kern %*% t(kern)
      kern_svd <- svd(kern_outer)
      kern_svd$d[-1] <- 0.0
      
      # reference standard deviation
      R2 <- R + diag(1e-8, nrow(R))
      S  <- solve(R2) 
      S2 <- S - S %*% kern %*% solve(t(kern) %*% S %*% kern) %*% t(kern) %*% S
      sigma_ref <- exp(sum(0.5 * log(diag(S2)))/ntimes)
      
      R_norm       <- R * sigma_ref^2 # precision matrix, normalized to have generalized marginal variance equal to 1
      R_norm_svd   <- svd(R_norm)
      R_norm_svd$d[(ntimes-1):ntimes] <- 0.0
      
      return(
            list(
                  D = D,
                  R = R,
                  R_svd = R_svd,
                  R_norm = R_norm,
                  R_norm_svd = R_norm_svd,
                  kern = kern,
                  kern_outer = kern_outer,
                  kern_svd = kern_svd,
                  sigma_ref = sigma_ref
            )
      )
}