<!-- README.md is generated from README.Rmd. Please edit that file -->
stemr - Stochastic Epidemic Models in R via Bayesian Data Augmentation
==========================================================================

[![Build Status](https://travis-ci.org/fintzij/stemr.svg?branch=master)](https://travis-ci.org/fintzij/stemr)

Fit stochastic epidemic models using Bayesian data augmentation. Currently implements approximate inference via the linear noise approximation with a suite of adaptive and fixed kernel MCMC algorithms for parameter proposals. Planned development includes an implementation of an agent-based data-driven data augmentation algorithm for efficient subject-level path proposals based on noisy aggregate counts at discrete observation times.
