# prior_from_quantiles
Computes mean and variance of a prior distribution from a quantile (as e.g. provided in Schmitt-Grohé/Uribe (2012): What's news in business cycles, Econometrica).

Run `run_prior_parameter_calculation.m` in Matlab after doing the necessary adjustments to the `prior_specification` cell. The file will generate a `prior_specification_Dynare.txt`file containing the prior specification using the mean and standard deviations as in Dynare.
