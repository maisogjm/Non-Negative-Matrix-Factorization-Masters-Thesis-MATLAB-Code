function muHat12sigma2Hat = gaussianParameterEstimatesMLE( vec , q )
% Input: vec - vector of eigenvalues (or eigenvalue equivalents)
%
%        q   - current 'guess' for number of components, an integer
%
% Output: muHat1    - estimate for Gaussian parameter mu1
%
%         muHat2    - estimate for Gaussian parameter mu2
%
%         sigma2Hat - estimate for Gaussian parameter sigma squared

% If desired, basic checks could be done here.  For example, q should be
% non-negative and less than the length of vec.  Since I know that the
% calling function will behave properly, I will dispense with the niceties.

% Compute MLEs of means.
vec1   = vec(1:q);
vec2   = vec(q+1:end);
muHat1 = mean(vec1);
if ( ~ isempty(vec2) )
    muHat2 = mean(vec2);
else
    muHat2 = 0;
end

% Compute SAMPLE variances (VAR(X,1) normalizes by N).
s2_1 = var(vec1,1);
if ( ~ isempty(vec2) )
    s2_2 = var(vec2,1);
else
    s2_2 = 0;
end

% Compute MLE of sigma squared, pooled estimate.
p         = length(vec);
sigma2Hat = ( ( (q-1) * s2_1 ) + ( (p-q-1) * s2_2 ) ) / (p-2);

muHat12sigma2Hat = [ muHat1 muHat2 sigma2Hat ];
return;
