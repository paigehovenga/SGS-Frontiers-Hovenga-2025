function [ar1_synth_nreal] = ar1_integration(x, nreal, lambda);

% Description:
% Creates nrealizations of AR1 process
% Used for computing confidence intervals and median values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% x = timeseries data
% nreal = number of bootstrapped realizations 
% lambda = 1/decorrelation timescale

% Outputs:
% ar1_synth_nreal = matrix of the bootstrapped AR1 synthetic time series

% Provenance: 
% Affiliation: Physical Sciences Laboratory, NOAA and Cooperative Institue for Research in
% Environmental Sciences, CU Boulder
% Based on: Sardeshmukh and Sura 2009 J. of Atmos. Sci.
% Written: Apr 2011 written by Gil Compo
% Adapted: to Matlab: Feb 2017 by Matt Newman
% Modified: Jun 2025 by Paige Hovenga (contact paige.hovenga@noaa.gov)

% Dependencies: 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b3ar1 = sqrt(var(x)*2*lambda);

dt = 0.01;%dt=.01;
len = length(x);
totalLen = round(len/dt);
ar1_synth_nreal = zeros(nreal,len);

% loop through number of bootstrap realizations 
for imc=1:nreal

    % create noise
    eta = randn(totalLen,2);
    eta1 = eta(:,1);

    xi(1) = 0;
    j=0;
    for i=1:totalLen;

        xi(i+1) = xi(i) - lambda*xi(i)*dt + b3ar1*eta(i)*sqrt(dt);

        if mod(i,1./dt)==0; % every .1 time unit, which is about 1 month given
            j=j+1;
            z3ar1save_tmp(j)=(xi(i)+xi(i+1))/2; %ar1_old
        end

    end

    ar1_synth_nreal(imc,:) = z3ar1save_tmp;

end

end