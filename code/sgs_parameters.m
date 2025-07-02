function [lambda, E, g, b, K, S] = sgs_parameters(x, lambda)

% Description:
% Fits the timeseries x to SGS model and solves for SGS parameters E, g, and b
% Reference: Sardeshmukh, P. D., G. P. Compo, and C. Penland, 2015:
% Need for Caution in Interpreting Extreme Weather Statistics.
% J. Climate, 28, 9166–9187, https://doi.org/10.1175/JCLI-D-15-0020.1.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% x = timeseries of data to compute the sgs for
% lambda = 1/decorrelation timescale (this is computed from 'compute_lambda_sgs.m' to deal with gaps in the data)

% Outputs:
% lambda = 1/decorrelation timescale (this is the same as the input)
% E, g and b = SGS parameters computed from the data
% K = excess kurtosis of the observations
% S = skewness of the observations

% Provenance: 
% Affiliation: Physical Sciences Laboratory, NOAA and Cooperative Institue for Research in
% Environmental Sciences, CU Boulder
% Based on: Sardeshmukh, P. D., G. P. Compo, and C. Penland, 2015
% Written: 
% Adapted: to Matlab: Feb 2017 by Matt Newman
% Modified: Jun 2025 by Paige Hovenga and Tongtong Xu (contact paige.hovenga@noaa.gov)

% Dependencies:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% small increment
kstep=.01;

% Based on Eqn (2)
len = length(x);
C = x*x'/(len-1);
C2 = x*(x.^2)'/(len-1);
C3 = x*(x.^3)'/(len-1);

V = sqrt(C); % sigma
S = C2/V^3; % skewness
K = C3/V^4 - 3; % excess kurtosis
eps = sqrt(1+(S^2/4));

% Check 'necessary and sufficient condition for both E2 and b2 to be positive'
b2=-1;
E2=-1;

while E2 < (eps-1)/eps || b2 < 0

    % Based on Eqn (8)
    E2 = 2/3*(K-(3/2)*S^2)/(K-S^2+2);
    g = S*V*(1-E2)/(2*sqrt(E2));
    b2 = 2*V^2*(1-E2/2-(1-E2)^2*S^2/(8*E2));

    K = K+kstep; %increase K if E2 and b2 do not meet conditions

end 

% result
E = sqrt(E2);
b = sqrt(b2);

end