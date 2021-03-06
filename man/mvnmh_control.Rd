% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mvnmh_control.R
\name{mvnmh_control}
\alias{mvnmh_control}
\title{Generate a list of settings for Metropolis-Hastings updates and adaptation
via the robust adaptive Metropolis algorithm (Vilhola, 2012).}
\usage{
mvnmh_control(
  n_updates = 1,
  target_acceptance = 0.234,
  max_scaling = Inf,
  scale_constant = 1,
  scale_cooling = 2/3,
  step_size = 1,
  stop_adaptation = 0,
  adaptation_offset = 0,
  nugget = 1e-05,
  nugget_cooling = 2/3,
  nugget_step_size = 1
)
}
\arguments{
\item{n_updates}{number of updates per iteration.}

\item{target_acceptance}{target acceptance rate, defaults to 0.234}

\item{max_scaling}{maximum global scaling factor, defaults to Inf.}

\item{scale_constant}{constant multiple of the adaptations determined by
\code{scale_cooling}.}

\item{scale_cooling}{rate at which to cool the adaptation, defaults to 0.5.
Adaptation contributions are governed by a harmonic sequence:
scale_constant/(iteration/step_size+1)^scale_cooling. The
\code{plot_adaptations} function may be used to plot the adaptation factors.}

\item{step_size}{adaptation increment for each iteration, defaults to 1.}

\item{stop_adaptation}{iteration at which adaptation should be terminated, 
defaults to 0 for no adaptation.}

\item{adaptation_offset}{iteration offset}

\item{nugget}{nugget for proposal covariance}

\item{nugget_cooling}{rate at which to cool the nugget, defaults to 0.9.
Adaptation contributions are governed by a harmonic sequence:
scale_constant/(iteration/step_size+1)^scale_cooling. The
\code{plot_adaptations} function may be used to plot the adaptation
factors.}

\item{nugget_step_size}{increment for each iteration, defaults to 100/number
of iterations.}
}
\value{
list with control settings for Metropolis-Hastings updates
}
\description{
Generate a list of settings for Metropolis-Hastings updates and adaptation
via the robust adaptive Metropolis algorithm (Vilhola, 2012).
}
