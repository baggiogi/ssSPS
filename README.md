# ssSPS
Matlab code to compute state-space SPS confidence regions of the paper

> G. Baggio, A. Carè, G. Pillonetto, "Finite-Sample Guarantees for State-Space System Identification", submitted, 2021.

The repository contains the script:

`main.m`: computes empirical coverage probability, component-wise box limits, mean distance from centroid, L2 ball radius of ssSPS regions and comparison with Gaussian Fisher (GF) and Observed Fisher (OF) regions for three noise scenarios described in the above paper.

and the following auxiliary functions:

`coverage.m`: computes empirical coverage probability of SPS, GF, OF regions.

`metrics.m`: computes performance metrics (component-wise box limits, mean distance from centroid, L2 ball radius) of SPS, GF, OF regions.

`vec.m`: vectorizes a matrix.
