MATLAB code for nonparametric Bayesian inference on the initial condition of the heat equation.

Author: Matteo Giordano, https://matteogiordano.weebly.com.

This repository is associated with the article "Bayesian inference for initial heat states with Gaussian series priors" by Matteo Giordano. The abstract of the paper is:

"We consider the statistical linear inverse problem of recovering the unknown initial heat state from noisy interior measurements over an inhomogeneous domain of the solution to the heat equation at a fixed time instant. We employ nonparametric Bayesian procedures with Gaussian series priors defined on the Dirichlet-Laplacian eigenbasis, yielding convenient conjugate posterior distributions with explicit expressions for posterior inference. We review recent theoretical results that provide asymptotic performance guarantees (in the large sample size limit) for the resulting posterior-based point estimation and uncertainty quantification. We further provide an implementation of the approach, and illustrate it via a numerical simulation study."

This repository contains the MATLAB code to replicate the numerical simulation study presented in Section 3 of the article. It contains the following:

1. GenerateObservations.m code to generate the observations (discrete point evaluations of the heat equation solution at a fixed time instant corrupted by additive Gaussian measurement errors).
2. pCNSeries.m code to implement posterior inference based on Gaussian series priors defined on the Dirichlet-Laplacian eigenbasis, via the Gaussian conjugate formulae.

For questions or for reporting bugs, please e-mail Matteo Giordano (matteo.giordano@unito.it).

Please cite the following article if you use this repository in your research: Giordano, M (2025). Bayesian inference for initial heat states with Gaussian series priors.
