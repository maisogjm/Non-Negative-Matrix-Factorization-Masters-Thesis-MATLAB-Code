function kHatvec = FYandBICMultiresStrategyKL( V , subsamplingRate )
%
% This function uses Fogel and Young's method and three BIC methods to estimate the number of components 'k'.
% Note that it calls Zhu and Ghodsi's method for determining the "elbow" in a scree plot.
% Also, it calls the MIT code to compute the NMF -- i.e., the Kullback-Liebler criterion.
%
% INPUTS:
% V              - data matrix.  Rows are observations, columns are variables.
% subsamplingRate - multiresolution factor for subsampling columns of V.  Must be at least 4; suggest 10.
% 
% OUTPUT:
%
% kHat    - estimate of the number of components 'k'
%

% Compute determinants
% Notations
% n: number of rows in original matrix
% p: number of column in original matrix
% np = n*p
% nc: number of NMF components
% W: Left factoring vectors
% H: Right factoring vectors
% S: Scaling factors
% X[0,k]: NMF Approximation of the original matrix on
% component k, reshaped into a column vector, and normalized.

% Make sure subsamplingRate is at least 4.
if ( subsamplingRate < 4 )
    disp('subsamplingRate must be at least 4!')
    kHat = -1;
    return
end

% Get dimensions of V.  Set the maximum possible number of components
% 'maxNC' equal to the minimum of [ numRows numCols ]; intuitively,
% maxNC should be at most rank(V).
[ n p ] = size(V);
maxNC   = min([ n p ] );

% Initialize Volume, BIC1, BIC2, BIC3, and alreadyTried vectors.
% Compute cmnFactors for the BIC measures.
% See Art Owen's paper, 'Bi-cross-validation of the SVD and the
% non-negative matrix factorization', page 20.
Volume       = zeros(maxNC,1); % Iterations over nc start at nc=2, so ignore first element.
BIC1         = zeros(maxNC,1); % Iterations over nc start at nc=2, so ignore first element.
BIC2         = zeros(maxNC,1); % Iterations over nc start at nc=2, so ignore first element.
BIC3         = zeros(maxNC,1); % Iterations over nc start at nc=2, so ignore first element.
rrssq        = zeros(maxNC,1); % Iterations over nc start at nc=2, so ignore first element.
alreadyTried = zeros(maxNC,1); % Set to 1 if Volume was computed for the corresponding index.
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

% Parameters for NTU NMF code.
tol       = 1e-3;
timelimit = 5*24*60*60; % Allow up to five days (?!)
maxiter   = 128000;

% Initialize bracketing points 'lowVal' and 'hiVal'.  Determine starting stepSize.
loVal_FY      = 2;
hiVal_FY      = maxNC;
stepSize_FY   = (hiVal_FY-loVal_FY)/subsamplingRate;

loVal_BIC1    = 2;
hiVal_BIC1    = maxNC;
stepSize_BIC1 = (hiVal_BIC1-loVal_BIC1)/subsamplingRate;

loVal_BIC2    = 2;
hiVal_BIC2    = maxNC;
stepSize_BIC2 = (hiVal_BIC2-loVal_BIC2)/subsamplingRate;

loVal_BIC3    = 2;
hiVal_BIC3    = maxNC;
stepSize_BIC3 = (hiVal_BIC3-loVal_BIC3)/subsamplingRate;

loVal_rrssq    = 2;
hiVal_rrssq    = maxNC;
stepSize_rrssq = (hiVal_rrssq-loVal_rrssq)/subsamplingRate;

% Loop over possible values for number of components.
% Start indexing nc at '2' instead of '1'; otherwise nmf_taiwan code misbehaves.
% Note: now using nmf_mit instead of nmf_taiwan.
np       = n*p;
finished = 0;

while ( finished == 0 )
    %-----------------------------------------------------------------------
    % Sampling based on Fogel and Young's method.
    %-----------------------------------------------------------------------

    % Determine the indices for which FY and BICs will be computed this round.
    % Subsample the interval [lowVal_FY, hiVal_FY] by a factor 'stepSize_FY'.
    indexList = round([0:subsamplingRate]*stepSize_FY)+loVal_FY;

    % Loop over the indices just computed, compute their Volumes.
    for nc=indexList
        if ( nc > length(alreadyTried) )
            continue
        end

        % If this value for 'k' hasn't been tried yet, evaluate it.
        if ( alreadyTried(nc) == 0 )
            % Scale factors -- set all to '1' for now.
            S = ones(nc,1);

            % Compute NMF.
            % [ W H ] = nmf_taiwan(V,rand(n,nc),rand(nc,p),tol,timelimit,maxiter);
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
            Volume(nc) = det(X(:,1:nc)'*X(:,1:nc));        % Volume[nc]=Det(X[0,1::nc]`*X[0,1::nc]);

            % Compute BIC1, BIC2, BIC3, and RRSSQ.
            Vhat      = W * H;
%            normDiff  = norm(Vhat-V,'fro');
%            logNorm   = log(normDiff*normDiff);
            logNorm   = RSSQ(V,Vhat);
            BIC1(nc)  = logNorm+(nc*cmnFactor1);
            BIC2(nc)  = logNorm+(nc*cmnFactor2);
            BIC3(nc)  = logNorm+(nc*cmnFactor3);
%            rrssq(nc) = RRSSQ(V,Vhat);
            rrssq(nc) = sqrt(logNorm); % Equivalent to RRSSQ(V,Vhat); saves recomputations.

            % Update alreadyTried(nc).
            alreadyTried(nc) = 1;
        end; % if alreadyTried
    end; % for nc

    %-----------------------------------------------------------------------
    % Sampling based on BIC1 method.
    %-----------------------------------------------------------------------

    % Determine the indices for which FY and BICs will be computed this round.
    % Subsample the interval [lowVal_BIC1, hiVal_BIC1] by a factor 'stepSize_BIC1'.
    indexList = round([0:subsamplingRate]*stepSize_BIC1)+loVal_BIC1;

    % Loop over the indices just computed, compute BIC1 for each.
    for nc=indexList
        if ( nc > length(alreadyTried) )
            continue
        end

        % If this value for 'k' hasn't been tried yet, evaluate it.
        if ( alreadyTried(nc) == 0 )
            % Scale factors -- set all to '1' for now.
            S = ones(nc,1);

            % Compute NMF.
            % [ W H ] = nmf_taiwan(V,rand(n,nc),rand(nc,p),tol,timelimit,maxiter);
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
            Volume(nc) = det(X(:,1:nc)'*X(:,1:nc));        % Volume[nc]=Det(X[0,1::nc]`*X[0,1::nc]);

            % Compute BIC1, BIC2, and BIC3.
            Vhat      = W * H;
%            normDiff  = norm(Vhat-V,'fro');
%            logNorm   = log(normDiff*normDiff);
            logNorm   = RSSQ(V,Vhat);
            BIC1(nc)  = logNorm+(nc*cmnFactor1);
            BIC2(nc)  = logNorm+(nc*cmnFactor2);
            BIC3(nc)  = logNorm+(nc*cmnFactor3);
%            rrssq(nc) = RRSSQ(V,Vhat);
            rrssq(nc) = sqrt(logNorm); % Equivalent to RRSSQ(V,Vhat); saves recomputations.

            % Update alreadyTried(nc).
            alreadyTried(nc) = 1;
        end; % if alreadyTried
    end; % for nc

    %-----------------------------------------------------------------------
    % Sampling based on BIC2 method.
    %-----------------------------------------------------------------------

    % Determine the indices for which FY and BICs will be computed this round.
    % Subsample the interval [lowVal_BIC2, hiVal_BIC2] by a factor 'stepSize_BIC2'.
    indexList = round([0:subsamplingRate]*stepSize_BIC2)+loVal_BIC2;

    % Loop over the indices just computed, compute BIC2 for each.
    for nc=indexList
        if ( nc > length(alreadyTried) )
            continue
        end

        % If this value for 'k' hasn't been tried yet, evaluate it.
        if ( alreadyTried(nc) == 0 )
            % Scale factors -- set all to '1' for now.
            S = ones(nc,1);

            % Compute NMF.
            % [ W H ] = nmf_taiwan(V,rand(n,nc),rand(nc,p),tol,timelimit,maxiter);
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
            Volume(nc) = det(X(:,1:nc)'*X(:,1:nc));        % Volume[nc]=Det(X[0,1::nc]`*X[0,1::nc]);

            % Compute BIC1, BIC2, BIC3, and RRSSQ.
            Vhat      = W * H;
%            normDiff  = norm(Vhat-V,'fro');
%            logNorm   = log(normDiff*normDiff);
            logNorm   = RSSQ(V,Vhat);
            BIC1(nc)  = logNorm+(nc*cmnFactor1);
            BIC2(nc)  = logNorm+(nc*cmnFactor2);
            BIC3(nc)  = logNorm+(nc*cmnFactor3);
%            rrssq(nc) = RRSSQ(V,Vhat);
            rrssq(nc) = sqrt(logNorm); % Equivalent to RRSSQ(V,Vhat); saves recomputations.

            % Update alreadyTried(nc).
            alreadyTried(nc) = 1;
        end; % if alreadyTried
    end; % for nc

    %-----------------------------------------------------------------------
    % Sampling based on BIC3 method.
    %-----------------------------------------------------------------------

    % Determine the indices for which FY and BICs will be computed this round.
    % Subsample the interval [lowVal_BIC3, hiVal_BIC3] by a factor 'stepSize_BIC3'.
    indexList = round([0:subsamplingRate]*stepSize_BIC3)+loVal_BIC3;

    % Loop over the indices just computed, compute BIC1 for each.
    for nc=indexList
        if ( nc > length(alreadyTried) )
            continue
        end

        % If this value for 'k' hasn't been tried yet, evaluate it.
        if ( alreadyTried(nc) == 0 )
            % Scale factors -- set all to '1' for now.
            S = ones(nc,1);

            % Compute NMF.
            % [ W H ] = nmf_taiwan(V,rand(n,nc),rand(nc,p),tol,timelimit,maxiter);
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
            Volume(nc) = det(X(:,1:nc)'*X(:,1:nc));        % Volume[nc]=Det(X[0,1::nc]`*X[0,1::nc]);

            % Compute BIC1, BIC2, BIC3, and RRSSQ.
            Vhat      = W * H;
%            normDiff  = norm(Vhat-V,'fro');
%            logNorm   = log(normDiff*normDiff);
            logNorm   = RSSQ(V,Vhat);
            BIC1(nc)  = logNorm+(nc*cmnFactor1);
            BIC2(nc)  = logNorm+(nc*cmnFactor2);
            BIC3(nc)  = logNorm+(nc*cmnFactor3);
%            rrssq(nc) = RRSSQ(V,Vhat);
            rrssq(nc) = sqrt(logNorm); % Equivalent to RRSSQ(V,Vhat); saves recomputations.

            % Update alreadyTried(nc).
            alreadyTried(nc) = 1;
        end; % if alreadyTried
    end; % for nc

    %-----------------------------------------------------------------------
    % Sampling based on Shao's Relative Root of Sum of Square Differences (RRSSQ) method.
    %-----------------------------------------------------------------------

    % Determine the indices for which FY and BICs will be computed this round.
    % Subsample the interval [lowVal_rrssq, hiVal_rrssq] by a factor 'stepSize_rrssq'.
    indexList = round([0:subsamplingRate]*stepSize_rrssq)+loVal_rrssq;

    % Loop over the indices just computed, compute BIC1 for each.
    for nc=indexList
        if ( nc > length(alreadyTried) )
            continue
        end

        % If this value for 'k' hasn't been tried yet, evaluate it.
        if ( alreadyTried(nc) == 0 )
            % Scale factors -- set all to '1' for now.
            S = ones(nc,1);

            % Compute NMF.
            % [ W H ] = nmf_taiwan(V,rand(n,nc),rand(nc,p),tol,timelimit,maxiter);
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
            Volume(nc) = det(X(:,1:nc)'*X(:,1:nc));        % Volume[nc]=Det(X[0,1::nc]`*X[0,1::nc]);

            % Compute BIC1, BIC2, BIC3, and RRSSQ.
            Vhat      = W * H;
%            normDiff  = norm(Vhat-V,'fro');
%            logNorm   = log(normDiff*normDiff);
            logNorm   = RSSQ(V,Vhat);
            BIC1(nc)  = logNorm+(nc*cmnFactor1);
            BIC2(nc)  = logNorm+(nc*cmnFactor2);
            BIC3(nc)  = logNorm+(nc*cmnFactor3);
%            rrssq(nc) = RRSSQ(V,Vhat);
            rrssq(nc) = sqrt(logNorm); % Equivalent to RRSSQ(V,Vhat); saves recomputations.

            % Update alreadyTried(nc).
            alreadyTried(nc) = 1;
        end; % if alreadyTried
    end; % for nc

    % Construct vector of all Volumes computed thus far, determine the Volume inflection point.
    alreadyTriedIndices = find(alreadyTried==1);
    volumeList = Volume(alreadyTriedIndices);

    % For Fogel and Young's method, use Zhu and Godhsi's method to update the current estimate for k.
    indexOfIndex_FY = zhuEstimate(volumeList);
    kHat_FY      = alreadyTriedIndices(indexOfIndex_FY);

    % Update the BIC estimates for k.
    [ Y1 indexOfIndex_BIC1 ] = min(BIC1(alreadyTriedIndices));
    kHat_BIC1             = alreadyTriedIndices(indexOfIndex_BIC1);

    [ Y2 indexOfIndex_BIC2 ] = min(BIC2(alreadyTriedIndices));
    kHat_BIC2             = alreadyTriedIndices(indexOfIndex_BIC2);

    [ Y3 indexOfIndex_BIC3 ] = min(BIC3(alreadyTriedIndices));
    kHat_BIC3             = alreadyTriedIndices(indexOfIndex_BIC3);

    [ Y4 indexOfIndex_rrssq ] = min(rrssq(alreadyTriedIndices));
    kHat_rrssq            = alreadyTriedIndices(indexOfIndex_rrssq);

    % If we have reached the finest stepsize, we are finished.
    if ( ( stepSize_FY    <= 1 ) ...
       & ( stepSize_BIC1  <= 1 ) ...
       & ( stepSize_BIC2  <= 1 ) ...
       & ( stepSize_BIC3  <= 1 ) ...
       & ( stepSize_rrssq <= 1 ) )
        finished = 1;

    % Else, proceed to the next finest resolution in the WHILE loop by reducing the stepSize
    % and narrowing the search range [loVal,hiVal].
    else
        numAlreadyTried = length(alreadyTriedIndices);

        %-------------------------------------------------------------------
        % Sampling based on Fogel and Young's method.
        %-------------------------------------------------------------------
        % With the current estimate for k (represented by indexOfIndex_FY), update loVal_FY.
        loIndexOfIndex = indexOfIndex_FY - 1;
        if ( loIndexOfIndex < 1 )
            loIndexOfIndex = 1;
        end
        loVal_FY = alreadyTriedIndices(loIndexOfIndex);

        % With the current estimate for k (represented by indexOfIndex_FY), update hiVal_FY.
        hiIndexOfIndex = indexOfIndex_FY + 1;
        if ( hiIndexOfIndex > numAlreadyTried )
            hiIndexOfIndex = numAlreadyTried;
        end
        hiVal_FY = alreadyTriedIndices(hiIndexOfIndex);

        % Update stepSize.
        stepSize_FY = (hiVal_FY-loVal_FY)/subsamplingRate;

        %-------------------------------------------------------------------
        % Sampling based on BIC1 method.
        %-------------------------------------------------------------------
        % With the current estimate for k (represented by indexOfIndex_BIC1), update loVal_BIC1.
        loIndexOfIndex = indexOfIndex_BIC1 - 1;
        if ( loIndexOfIndex < 1 )
            loIndexOfIndex = 1;
        end
        loVal_BIC1 = alreadyTriedIndices(loIndexOfIndex);

        % With the current estimate for k (represented by indexOfIndex_BIC1), update hiVal_BIC1.
        hiIndexOfIndex = indexOfIndex_BIC1 + 1;
        if ( hiIndexOfIndex > numAlreadyTried )
            hiIndexOfIndex = numAlreadyTried;
        end
        hiVal_BIC1 = alreadyTriedIndices(hiIndexOfIndex);

        % Update stepSize.
        stepSize_BIC1 = (hiVal_BIC1-loVal_BIC1)/subsamplingRate;

        %-------------------------------------------------------------------
        % Sampling based on BIC2 method.
        %-------------------------------------------------------------------
        % With the current estimate for k (represented by indexOfIndex_BIC2), update loVal_BIC2.
        loIndexOfIndex = indexOfIndex_BIC2 - 1;
        if ( loIndexOfIndex < 1 )
            loIndexOfIndex = 1;
        end
        loVal_BIC2 = alreadyTriedIndices(loIndexOfIndex);

        % With the current estimate for k (represented by indexOfIndex_BIC2), update hiVal_BIC2.
        hiIndexOfIndex = indexOfIndex_BIC2 + 1;
        if ( hiIndexOfIndex > numAlreadyTried )
            hiIndexOfIndex = numAlreadyTried;
        end
        hiVal_BIC2 = alreadyTriedIndices(hiIndexOfIndex);

        % Update stepSize.
        stepSize_BIC2 = (hiVal_BIC2-loVal_BIC2)/subsamplingRate;

        %-------------------------------------------------------------------
        % Sampling based on BIC3 method.
        %-------------------------------------------------------------------
        % With the current estimate for k (represented by indexOfIndex_BIC3), update loVal_BIC3.
        loIndexOfIndex = indexOfIndex_BIC3 - 1;
        if ( loIndexOfIndex < 1 )
            loIndexOfIndex = 1;
        end
        loVal_BIC3 = alreadyTriedIndices(loIndexOfIndex);

        % With the current estimate for k (represented by indexOfIndex_BIC3), update hiVal_BIC3.
        hiIndexOfIndex = indexOfIndex_BIC3 + 1;
        if ( hiIndexOfIndex > numAlreadyTried )
            hiIndexOfIndex = numAlreadyTried;
        end
        hiVal_BIC3 = alreadyTriedIndices(hiIndexOfIndex);

        % Update stepSize.
        stepSize_BIC3 = (hiVal_BIC3-loVal_BIC3)/subsamplingRate;

        %-------------------------------------------------------------------
        % Sampling based on Shao's RRSSQ method.
        %-------------------------------------------------------------------
        % With the current estimate for k (represented by indexOfIndex_rrssq), update loVal_rrssq.
        loIndexOfIndex = indexOfIndex_rrssq - 1;
        if ( loIndexOfIndex < 1 )
            loIndexOfIndex = 1;
        end
        loVal_rrssq = alreadyTriedIndices(loIndexOfIndex);

        % With the current estimate for k (represented by indexOfIndex_rrssq), update hiVal_rrssq.
        hiIndexOfIndex = indexOfIndex_rrssq + 1;
        if ( hiIndexOfIndex > numAlreadyTried )
            hiIndexOfIndex = numAlreadyTried;
        end
        hiVal_rrssq = alreadyTriedIndices(hiIndexOfIndex);

        % Update stepSize.
        stepSize_rrssq = (hiVal_rrssq-loVal_rrssq)/subsamplingRate;
    end; % If stepSize_FY stepSize_BIC1 stepSize_BIC2 stepSize_BIC3 stepSize_rrssq
end; % while finished

% Return vector of estimated k's.
kHatvec = [ kHat_FY kHat_BIC1 kHat_BIC2 kHat_BIC3 kHat_rrssq ];

return
