
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stemr <img src="man/figures/logo.png" align="right" alt="" width="120" />

## Bayesian Inference for Stochastic Epidemic Models via the Linear Noise Approximation

The **stemr** R package implements a Bayesian data augmentation
algorithm for fitting stochastic epidemic models (SEMs) with arbitrary
dynamics to partially observed epidemic count data via the linear noise
approximation (LNA). This repository and package exist to expose and
archive the code used in the paper, “A Linear Noise Approximation for
Stochastic Epidemic Models Fit to Partially Observed Incidence Counts”
by J. Fintzi, J. Wakefield, and V.N. Minin. The package provides
facilities for exact simulation of SEMs via Gillespie’s direct
algorithm, and for simulation of approximate epidemic paths via the LNA.
Inference is carried out using a Bayesian data augmentation algorithm.
LNA paths are sampled via the elliptical slice sampling algorithm of
Murray et al. (2010), and the package offers several choices for
adaptive and fixed kernel MCMC algorithms for parameter updates. The
package also allows for time varying covariates, time varying
parameters, and deterministic forcings. These functionalities are not
yet fully documented, but sample code is available in the
‘/.development\_files’ directory on the GitHub repository.

## Package installation

To install the `stemr` package, clone this repository and build the
package from sources. You should be able to rebuild in the usual way
once you clone the package repo and install the other dependencies
(`odeintr`, `extraDistr`, `Rcpp`, `RcppArmadillo`, and `BH`).
Computationally intensive components in `stemr` package are implemented
in C++. Hence, it is important to also make sure that your C++ toolchain
is set up properly, e.g., by following instructions given in the
[Stan](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started)
documentation, and that Rtools has been added to your system path. If
you are working on a Windows machine, you may need to take additional
steps to ensure your toolset is in order.

## Vignettes and Simulation Code

The code for reproducing the simulations in Fintzi, et. al (2020) can be
found in the
[‘/.development\_files/incidence\_lna\_simulations’](https://github.com/fintzij/stemr/tree/master/.development_files/incidence_lna_simulations)
folder. There are two vignettes included in this package to help
familiarize users with its basic functionality and that reproduce models
the SEMs fit via the LNA in Fintzi, et al. (2020). The
[stemr](https://github.com/fintzij/stemr/blob/master/vignettes/stemr.Rmd)
vignette provides an introduction to the `stemr` package, and
demonstrates how to simulate from and fit an SIR model via the LNA and
ODE. The
[ebola\_westafrica](https://github.com/fintzij/stemr/blob/master/vignettes/ebola_westafrica.Rmd)
vignette demonstrates how to simulate from and fit a multi-country model
for Ebola transmission, and then provides code to fit the model to data
from the 2014-2015 outbreak in West Africa.
