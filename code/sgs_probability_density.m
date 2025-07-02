function [sgspdf] = sgs_probability_density (Ec,gc,bc,x)

% Description:
% Determines SGS ("Stochastically-Generated Skewed") distribution given input range of 
% values x (assumed standardized). Based on Sardeshmukh and Sura 2009 J. of Atmos. Sci.

% Assume that the correlated additive and multiplicative noise differential equation is:
% dx/dt = -(lambda + 0.5*Ec^2)x+ bc*eta1+(Ec*x+gc)*eta2 - 0.5*Ec*gc

% Probability density function of this equation is given in SS09 as:
% P(x) = 1/K * ((Ec*x+gc)^2  + bc^2)^(-(1+lambda/Ec^2))*exp(2*gc/bc*(lambda/Ec^2)*atan((Ec*x+gc)/bc))
% where K is the constant as derived by C. Penland:
% K = 2*!pi/Ec*(2*b)^(2*lambda/Ec^2+1)*% Gamma(2*lambda/Ec^2+1)/(Gamma(lambda/Ec^2+1 - i * gc/bc*(lambda/Ec^2))*Gamma(lambda/Ec^2+1 + i * gc/bc*(lambda/Ec^2)) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Inputs:
% Ec = SGS constant E
% gc = SGS constant g
% bc = SGS constant b
% x = vector of  x values to be evaluated

% Outputs:
% sgspdf = the probability density of the SGS PDF

% Other variables:
% eta1 = white noise of unit amplitude from an N(0,1) distribution
% eta2 = white noise of unit amplitude from an N(0,1) distribution
% lambda = decorrelation timescale 

% Provenance: 
% Affiliation: Physical Sciences Laboratory, NOAA and Cooperative Institue for Research in
% Environmental Sciences, CU Boulder
% Based on: Sardeshmukh and Sura 2009 J. of Atmos. Sci.
% Written: Apr 2011 written by Gil Compo
% Adapted: to Matlab: Feb 2017 by Matt Newman
% Modified: Jun 2025 by Paige Hovenga (contact paige.hovenga@noaa.gov)

% Dependencies: 
% gammalog.m 
% we evaluate the logarithm of the SGS PDF using the more precise Lngamma and then take the exponential 
% matlab built-in lngamma is real only, so here we use gammalog by Paul Godfrey (see matlab file exchange)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% lambda = decorrelation timescale 
lambda=1; % always, due to normalization

dpi=3.14159265358979323846264338327950288419; % so many digits!!

i=0;
for xi = x

  i=i+1;
  
  nu = lambda/(Ec^2.0d0);

  com_valp = complex(nu+1.0d0,gc/bc*nu);
  com_valm = complex(nu+1.0d0,-1.0d0*gc/bc*nu);
  
  term1 = log(Ec/(2*dpi));
  term2 = (2.0d0*nu+1.0d0)*log(2.0d0*bc);
  term3 = gammalog(com_valm);
  
  term4 = gammalog(com_valp);
  
  term5 = -1.0d0*gammalog(2.0d0*nu+1.0d0);
  xparent=(Ec*xi+gc)^2  + bc^2;
  term6 = -1.0d0*(1.0d0+nu)*log(xparent);
  term7 = 2*gc/bc*(nu)*atan((Ec*xi+gc)/bc);
  sum_arg = term1 + term2 + real(term3+term4)+term5+term6+term7;
  
  sgspdf(i) = exp(sum_arg);

end
