function [ Volume BIC1 BIC2 BIC3 rrssq ] =...
    FY_BIC_RRSSQ_KL(V,indexList)
% Get dimensions of V.  Set the maximum possible number of components
% 'maxNC' equal to the minimum of [ numRows numCols ]; intuitively,
% maxNC should be at most rank(V).
[ n p ] = size(V);
numK    = length(indexList);
np      = n*p;

% Initialize Volume, BIC1, BIC2, BIC3, and rrssq output vectors.
Volume       = zeros(numK,1); % Iterations over nc start at nc=2, so ignore first element.
BIC1         = zeros(numK,1); % Iterations over nc start at nc=2, so ignore first element.
BIC2         = zeros(numK,1); % Iterations over nc start at nc=2, so ignore first element.
BIC3         = zeros(numK,1); % Iterations over nc start at nc=2, so ignore first element.
rrssq        = zeros(numK,1); % Iterations over nc start at nc=2, so ignore first element.

% Compute cmnFactors for the BIC measures.
% See Art Owen's paper, 'Bi-cross-validation of the SVD and the
% non-negative matrix factorization', page 20.
sm = sqrt(n);
sn = sqrt(p);
if ( sm < sn )
    c = sm;
else
    c = sn;
end
mnFactor   = ( n + p ) / ( n * p );
cmnFactor1 = mnFactor * log(1/mnFactor);
cmnFactor2 = mnFactor*log(c*c);
cmnFactor3 = cmnFactor2 / (c*c);

% Loop over the indexList
outIndex = 1; % index into output vectors
for nc=indexList

    % Scale factors -- set all to '1' for now.
    S = ones(nc,1);

    % Compute NMF.
    [ W H ] = nmf_mit(V,nc,0);

    % Fogel and Young's method - MATLAB code.      % ORIGINAL JSL CODE BELOW
    X = zeros(np,nc);                              % X=J(np,nc,.); % Initialize approximation
    for k=1:nc                                     % for (k=1,k<=nc,k++,
        X(:,k) = reshape(S(k)*W(:,k)*H(k,:),np,1); %    X[0,k]=shape(S[k]*W[0,k]*H[0,k]`,np,1);
        denom  = sqrt(X(:,k)'*X(:,k));
        if ( denom ~= 0 )
            X(:,k) = X(:,k)/denom;                 %    X[0,k]/=sqrt(X[0,k]`*X[0,k]);
        end
    end                                            % );%  end for k
    Volume(outIndex) = det(X(:,1:nc)'*X(:,1:nc));     % Volume[nc]=Det(X[0,1::nc]`*X[0,1::nc]);

    % Compute BIC1, BIC2, BIC3, and RRSSQ.
    Vhat      = W * H;
    % normDiff  = norm(Vhat-V,'fro');
    % logNorm   = log(normDiff*normDiff);
    logNorm   = RSSQ(V,Vhat);
    BIC1(outIndex)  = logNorm+(nc*cmnFactor1);
    BIC2(outIndex)  = logNorm+(nc*cmnFactor2);
    BIC3(outIndex)  = logNorm+(nc*cmnFactor3);
    rrssq(outIndex) = sqrt(logNorm); % Equivalent to RRSSQ(V,Vhat); saves recomputations.
    % rrssq(nc) = RRSSQ(V,Vhat);

    % Increment output vector index.
    outIndex = outIndex + 1;
end; % for nc

