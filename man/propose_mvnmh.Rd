% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{propose_mvnmh}
\alias{propose_mvnmh}
\title{Multivariate normal Metropolis-Hastings proposal}
\usage{
propose_mvnmh(params_prop, params_cur, kernel_cov_chol, nugget)
}
\arguments{
\item{params_prop}{vector in which the proposed parameters should be stored}

\item{params_cur}{vector containing the current parameter vector}

\item{kernel_cov_chol}{cholesky of the kernel covariance}

\item{nugget}{zero if adaptation is not ongoing}
}
\value{
propose new parameter values in place
}
\description{
Multivariate normal Metropolis-Hastings proposal
}
