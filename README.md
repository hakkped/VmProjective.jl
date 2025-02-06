# VmProjective.jl

# VmProjective

<!-- [![Build Status](https://github.com/hakkped/FractalDim.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/hakkped/FractalDim.jl/actions/workflows/CI.yml?query=branch%3Amain) -->

Compute the [co-moving velocity](https://link.springer.com/article/10.1007/s11242-022-01783-7) from a projective constitutive relation. A link to an article using the routine will be added soon. The code might be turned into a proper package in the future.

The code serves as an example of the routines involved, however, they are quite simple. It uses functions from [ImageProjectiveGeometry.jl](https://github.com/peterkovesi/ImageProjectiveGeometry.jl) to compute a 1d homography via [Singular Value Decomposition (SVD)](https://en.wikipedia.org/wiki/Singular_value_decomposition).

Input data are vectors of saturation and the two fluid velocities, in addition to the values of some physical parameters depending on the input data. This can of course be changed to ones liking in the code. An example is provided in `runTestExample.jl` for the relperm-data. One simply specifies a file path, define the relevant physical quantities, fetch the data with `getrelpermdata` and compute results using `computeall`. Plotting examples are included.
