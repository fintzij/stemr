% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{pars2lnapars}
\alias{pars2lnapars}
\title{Insert parameters into each row of a parameter matrix}
\usage{
pars2lnapars(lnapars, parameters)
}
\arguments{
\item{lnapars}{matrix of lna parameters, constants, and time-varying covars}

\item{parameters}{vector of parameters to be copied into the matrix}
}
\value{
modifies the lna parameter matrix in place
}
\description{
Insert parameters into each row of a parameter matrix
}
