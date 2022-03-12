# State-space SPS (ssSPS)
Matlab code to compute state-space SPS confidence regions as described in the paper:

> G. Baggio, A. Carè, G. Pillonetto, "Finite-Sample Guarantees for State-Space System Identification", submitted, 2022.

The repository contains:

- the script `main.m`: computes empirical coverage probability, component-wise box limits, mean distance from centroid, L2 ball radius of ssSPS regions and comparison with Gaussian Fisher (GF) and Observed Fisher (OF) regions for three noise scenarios described in the above paper.

- the auxiliary function `coverage.m`: computes empirical coverage probability of SPS, GF, OF regions.

- the auxiliary function `metrics.m`: computes performance metrics (component-wise box limits, mean distance from centroid, L2 ball radius) of SPS, GF, OF regions.

- the auxiliary function`vec.m`: vectorizes a matrix.

***

Author : G. Baggio <br/>
E-mail : baggio [dot] giacomo [at] gmail [dot] com
