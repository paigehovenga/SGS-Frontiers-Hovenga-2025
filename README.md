# SGS-Frontiers-Hovenga-2025
Stochastically Generated Skewed (SGS) code and data files associated with Hovenga et al., 2025

Code (SGS and AR1)
This set of functions can be used to solve for SGS parameters, compute the associated PDF, generate synthetic data, and compute confidence intervals
The order to use them are:
  1. sgs_compute_lambda.m - this will compute the 1/decorrelation timescale of the timeseries observations
  2. sgs_parameters.m - computes the SGS coefficients E, g, and b
  3. sgs_probability_density.m - computes the SGS PDF
  4. sgs_integration.m - creates nrealizations of SGS synthetic data
  5. ar1_integration.m - creates nrealizations of the AR1 (Gaussian) data
  6. estimate_pdf_CI.m - estimates the median and confidence intervals from either the SGS or AR1 nrealizations

Data
This dataset includes the computed lambda and SGS parameters E, g, and b from the 148 tide gauges associated with Hovenga et al., 2015
