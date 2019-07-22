<!-- README.md is generated from README.Rmd. Please edit that file -->
stemr
=====

Baysian Inference for Stochastic Epidemic Models via the Linear Noise Approximation
-----------------------------------------------------------------------------------

The **stemr** R package implements a Bayesian data augmentation
algorithm for fitting stochastic epidemic models (SEMs) with arbitrary
dynamics to partially observed epidemic count data via the linear noise
approximation (LNA). This repository and package exist to expose and
archive the code used in the paper, “Outbreak Modeling via the Linear
Noise Approximation: Fast Approximations for Stochastic Epidemic Models
Fit to Partially Observed Incidence Counts” by J. Fintzi, J. Wakefield,
and V.N. Minin. The package provides facilities for exact simulation of
SEMs via Gillespie’s direct algorithm, and for simulation of approximate
epidemic paths via the LNA. Inference is carried out using a Bayesian
data augmentation algorithm. LNA paths are sampled via the elliptical
slice sampling algorithm of Murray et al. (2010), and the package offers
several choices for adaptive and fixed kernel MCMC algorithms for
parameter updates. The package also allows for time varying covariates,
time varying parameters, and deterministic forcings.

Getting started
---------------

This package may be installed directly from GitHub using the
**devtools** package. **IMPORTANT** It is critical that the package is
installed with byte compilation disabled. The package will not work
properly when byte compiled, which is the default since R 3.4.

    library(devtools)
    install_github("fintzij/stemr", build_vignettes=TRUE, build_opts = c("--no-byte-compile")) 
    library(stemr)

There are several vignettes included in this package. The first diagrams
the basic functionality of the package. The others contain walk-throughs
of the simulations and analyses presented in the paper. These are
accessible as follows:

    vignette("stemr")
    vignette("sir_coverage")
    vignette("ebola_westafrica")

We also provide vignettes demonstrating how to use the `pomp` package to
fit the models in the coverage simulation and analysis of data from the
2013-2015 outbreak of Ebola in West Africa. These may be accessed as
follows:

    vignette("sir_pomp")
    vignette("ebola_pomp")
