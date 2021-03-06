% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/forcing.R
\name{forcing}
\alias{forcing}
\title{Declare a time varying covariate to be a forcing variable that moves
individuals between model compartments at discrete times. Flow in and out of
model compartments is allocated proportionally to the compartment counts in
the source and destination compartments.}
\usage{
forcing(tcovar_name, from, to)
}
\arguments{
\item{tcovar_name}{name of the time varying covariate}

\item{from}{vector of compartment names from which individuals exit}

\item{to}{vector of compartment names to which individuals enter}
}
\value{
list containing the forcing specification
}
\description{
Declare a time varying covariate to be a forcing variable that moves
individuals between model compartments at discrete times. Flow in and out of
model compartments is allocated proportionally to the compartment counts in
the source and destination compartments.
}
