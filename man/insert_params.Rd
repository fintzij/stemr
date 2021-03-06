% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/insert_params.R
\name{insert_params}
\alias{insert_params}
\title{Insert natural scale parameters into a parameter matrix}
\usage{
insert_params(parmat, param_blocks, nat, prop, rowind)
}
\arguments{
\item{parmat}{matrix into which parameters should be inserted}

\item{param_blocks}{list of parameter blocks}

\item{nat}{should parameters be inserted on their estimation scale}

\item{prop}{should proposed parameters be inserted?}

\item{rowind}{C++ row index}
}
\value{
insert natural scale parameters into a parameter matrix
}
\description{
Insert natural scale parameters into a parameter matrix
}
