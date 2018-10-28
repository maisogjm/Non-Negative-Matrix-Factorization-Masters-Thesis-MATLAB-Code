function lq = logLikelihood( vec , q )
% Input: vec - vector of eigenvalues (or eigenvalue equivalents)
%
%        q   - current 'guess' for number of components, an integer
%
% Output: lq - log likelihood for q

% Obtain MLE estimates for Gaussian parameters.
%[ muHat1 muHat2 sigma2Hat ] = gaussianParameterEstimatesMLE(vec,q);
muHat12sigma2Hat = gaussianParameterEstimatesMLE(vec,q);
muHat1    = muHat12sigma2Hat(1);
muHat2    = muHat12sigma2Hat(2);
sigma2Hat = muHat12sigma2Hat(3);

% First summation term.
lq = 0;
for i=1:q
    lq = lq + log(f(vec(i),muHat1,sigma2Hat));
end

% Second summation term.
p = length(vec);
for j=q+1:p
    lq = lq + log(f(vec(j),muHat2,sigma2Hat));
end
return;
