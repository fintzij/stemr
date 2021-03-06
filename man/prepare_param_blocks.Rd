% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_param_blocks.R
\name{prepare_param_blocks}
\alias{prepare_param_blocks}
\title{Prepare the parameter blocks for MCMC}
\usage{
prepare_param_blocks(
  param_blocks,
  parameters,
  param_codes,
  iterations,
  mcmc_restart
)
}
\arguments{
\item{param_blocks}{list of parameter blocks}

\item{parameters}{vector of parameters}

\item{param_codes}{vector of parameter codes}

\item{iterations}{number of MCMC iterations}

\item{mcmc_restart}{is mcmc being restarted}
}
\value{
validated param_block list with mcmc bookkeeping objects
}
\description{
Prepare the parameter blocks for MCMC
}
