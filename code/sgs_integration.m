function sgs_synth_nreal = sgs_integration(x, nreal, lambda)

% Description:
% Creates nrealizations of SGS synthetic data 
% Used for computing confidence intervals and median values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% x = timeseries data
% nreal = number of bootstrapped realizations 
% lambda = 1/decorrelation timescale

% Outputs:
% sgs_synth_nreal = matrix of the bootstrapped SGS synthetic time series

% Provenance: 
% Affiliation: Physical Sciences Laboratory, NOAA and Cooperative Institue for Research in
% Environmental Sciences, CU Boulder
% Based on: Sardeshmukh and Sura 2009 J. of Atmos. Sci.
% Written: Apr 2011 written by Gil Compo
% Adapted: to Matlab: Feb 2017 by Matt Newman
% Modified: Jun 2025 by Paige Hovenga (contact paige.hovenga@noaa.gov)

% Dependencies: 
% sgs_parameters.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% parameters of sgs distribution
[lambda,E,g,b] = sgs_parameters(x, lambda);


dt=.01;
len = length(x);
totalLen = round(len/dt); 
sgs_synth_nreal = zeros(nreal,len);

% loop through number of bootstrap realizations 
for imc=1:nreal 
    
    % create noise
    eta = randn(totalLen,2);
    eta1 = eta(:,1);
    eta2 = eta(:,2);

    xi=0;
    j=0;
    for i=1:totalLen        
        
        prod1 = -(1 + 0.5*E^2)*xi;
        prod2 = b*eta1(i) + (E*xi+g)*eta2(i);
        prod3 = -0.5*E*g ; 		% drift
        xt = xi +  prod1*dt*lambda + prod2*sqrt(dt*lambda) + prod3*dt*lambda;     
        
        prod1 = -0.5*(1 + 0.5*E^2)*(xi+xt);
        prod2 = b*eta1(i) + (E*(xi+xt)*0.5+g)*eta2(i);
        prod3 = -0.5*E*g ; 		% drift
        xp = xi +  prod1*dt*lambda + prod2*sqrt(dt*lambda) + prod3*dt*lambda;
        
        xi = xp;
        if mod(i,1./dt)==0 % every .1 time unit, which is about 1 month given
            j = j + 1;
            z3(j) = xp;
        end
    end
    
    sgs_synth_nreal(imc,:) = z3;
end

end