function pointDensity = f( d , mu , sigma2 )
% Input: d      - realization from the Gaussian distribution
%
%        mu     - mean of the Gaussian distribution
%
%        sigma2 - variance of the Gaussian distribution
%
% Output: pointD - point density of di for parameters mu and sigma2

delta        = d - mu;
delta2       = delta*delta;
twiceSigma2  = 2*sigma2;
pointDensity = (1/sqrt(twiceSigma2*pi)) * exp( (-delta2) / twiceSigma2 );
return;
