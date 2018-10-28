function kHat = fogelAndYoungMultiresStrategy( V , subsamplingRate )
%
% This function uses Fogel and Young's method to estimate the number of components 'k'.
% Note that it calls Zhu and Ghodsi's method for determining the "elbow" in a scree plot.
% Local variables 'a', 'b', and 'c' are a "bracketing triplet"; see Section 10.1
% of Numerical Recipes in C (2nd ed., 1992) for the idea.
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
maxNC = min([ n p ] );

% Initialize Volume, kVector vectors, index vectors, volume.
Volume       = zeros(maxNC,1); % Iterations over nc start at nc=2, so ignore first element.
kVector      = zeros(maxNC,1); % Estimates of k     start at nc=4, so ignore first three elements.
alreadyTried = zeros(maxNC,1); % Set to 1 if Volume was computed for the corresponding index.
baseVec      = 0:subsamplingRate;

% Parameters for NTU NMF code.
tol       = 1e-3;
timelimit = 5*24*60*60; % Allow up to five days (?!)
maxiter   = 128000;

% Initialize bracketing points 'lowVal' and 'hiVal'.  Determine starting stepSize.
loVal    = 2;
hiVal    = maxNC;
stepSize = (hiVal-loVal)/subsamplingRate;

% Loop over possible values for number of components.
% Start indexing nc at '2' instead of '1'; otherwise nmf_taiwan code misbehaves.
np       = n*p;
finished = 0;

while ( finished == 0 )
    % Determine the indices for which Volume will be computed this round.
    % Subsample the interval [lowVal, hiVal] by a factor 'stepSize'.
    indexList = round([0:subsamplingRate]*stepSize)+loVal;

    % Loop over the indices just computed, compute their Volumes.
    for nc=indexList

        % If this value for 'k' hasn't been tried yet, evaluate it.
        if ( alreadyTried(nc) == 0 )
            % Scale factors -- set all to '1' for now.
            S = ones(nc,1);

            % Compute NMF.
            [ W H ] = nmf_taiwan(V,rand(n,nc),rand(nc,p),tol,timelimit,maxiter);

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

            % Update alreadyTried(nc).
            alreadyTried(nc) = 1;
        end; % if alreadyTried
    end; % for nc

    % Construct vector of all Volumes computed thus far, determine the Volume inflection point.
    alreadyTriedIndices = find(alreadyTried==1);
    volumeList = Volume(alreadyTriedIndices);

    % Use Zhu and Godhsi's method to update the current estimate for k.
    indexOfIndex = zhuEstimate(volumeList);
    newKhat      = alreadyTriedIndices(indexOfIndex);

    % If we have reached the finest stepsize, we are finished.
    if ( stepSize <= 1 )
        finished = 1;

    % Else, proceed to the next finest resolution in the WHILE loop by reducing the stepSize.
    else
        numAlreadyTried = length(alreadyTriedIndices);

        % With the current estimate for k (represented by indexOfIndex), update loVal.
        loIndexOfIndex = indexOfIndex - 1;
        if ( loIndexOfIndex < 1 )
            loIndexOfIndex = 1;
        end
        loVal = alreadyTriedIndices(loIndexOfIndex);

        % With the current estimate for k (represented by indexOfIndex), update hiVal.
        hiIndexOfIndex = indexOfIndex + 1;
        if ( hiIndexOfIndex > numAlreadyTried )
            hiIndexOfIndex = numAlreadyTried;
        end
        hiVal = alreadyTriedIndices(hiIndexOfIndex);

        % Update stepSize.
        stepSize = (hiVal-loVal)/subsamplingRate;
    end; % If stepSize
end; % while finished

% Return last estimate of k.
kHat = newKhat;

%disp(alreadyTriedIndices)
%plot(alreadyTriedIndices,Volume(alreadyTriedIndices),'*-')

return
