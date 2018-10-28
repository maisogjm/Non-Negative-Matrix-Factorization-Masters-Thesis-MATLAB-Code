function [ qHat ] = zhuEstimate( vec )
% Reference: Zhu M and Ghodsi A, Comput Statist (2006); 51:817-930.
%
% Input: vec - vector of eigenvalues (or eigenvalue equivalents)
%
% Output: qHat - estimate for the number of components

% tmpQ is only for plotting, remove for final code
tmpQ = zeros(length(vec),1);

% Make sure that vec is sorted
% vec = sort(vec(:),1,'descend');
vec  = sort(vec(:),1);
nVec = length(vec);
vec  = vec(nVec:-1:1);

% Find 'q' where the log likelihood is maximal.
p          = length(vec);
qHat       = NaN;  %Initialize
maxLogLike = -Inf; %Initialize
for q=1:p
    tmpQ(q) = logLikelihood(vec,q);
    if ( tmpQ(q) > maxLogLike )
        maxLogLike = tmpQ(q);
        qHat = q; % Keep track of maxarg.
    end
end

% save /tmp/vec.mat vec tmpQ

% Comment out the following plot commands when this code is finalized.
% if ( 0 )
% subplot(1,2,1)
% bar(vec)
% title('Scree Plot')
% xlabel('Dimension')
% subplot(1,2,2)
% plot(tmpQ)
% title('Profile Log-Likelihood')
% xlabel('Dimension')
% end

return;
