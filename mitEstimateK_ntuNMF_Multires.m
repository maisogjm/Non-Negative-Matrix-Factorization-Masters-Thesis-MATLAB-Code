function [ kHat, cophCoefVector, alreadyTried ] = mitEstimateK_ntuNMF_Multires( V , subsamplingRate , cophThreshFactor , numRuns )
%
% This function uses Brunet et al.'s method to estimate the number of components, using a subsampling
% strategy in an attempt to optimize kHat without doing an exhaustive search.
% Note that we are using the National Taiwan University NMF code to compute the NMF, because it's fast.
%
% INPUTS:
% V               - data matrix.  Rows are observations, columns are variables.
% subsamplingRate - multiresolution factor for subsampling columns of V.  Must be at least 4; suggest 10.
% cophThreshFactor   - factor for computing the minimum threshold for the cophenetic coefficient.  The smallest 'kHat' such that the
%                   cophentic coefficient is greater than this threshold will be returned as the
%                   estimate for 'k'.  The proper choice will vary from data set to data set!
%                   Note: sometimes all of the cophenetic coefficients are
%                   subthreshold!  To address this problem, the threshold will
%                   be computed by multiplying cophThreshFactor by the maximum
%                   cophenetic correlation computed thus far.
% numRuns         - number of runs to average over for signal averaging, passed as the 'nloop' argument
%                   to nmfconsensus().  Obviously, the bigger this number the more stable and
%                   reproducible the estimate, but it will take longer to compute.
% 
% OUTPUT:
%
% kHat           - estimate of the number of components 'k'
% cophCoefVector - vector of cophenetic coefficients
% alreadyTried   - vector indicating which values of kHat were actually assessed
%

VERBOSE = 0;

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

% Initialize cophCoefVector and alreadyTried.
cophCoefVector = zeros(maxNC,1); % Iterations over kHatc start at kHatc=2, so ignore first element.
alreadyTried   = zeros(maxNC,1); % Set to 1 if cophenetic coefficient was computed for the corresponding index.

% Evaluate cophenetic coefficient for kHat = 2 thru 10, just to "catch" a few high values for starters.
% Note transposition of the data matrix, because Brunet's convention is that
% rows are variables and columns are observations.
% Don't start at kstart=1 because the NTU NMF code doesn't like it.
kstart = 2;
kend   = 6;

if ( VERBOSE )
    disp(sprintf('mitEstimateK_ntuNMF_Multires: kstart = %d\n',kstart))
    disp(sprintf('mitEstimateK_ntuNMF_Multires: kend = %d\n',kend))
end

consensus                            = nmfconsensus_taiwan(V',kstart,kend,numRuns,0);
[ ordcons, clustid, ordindex, coph ] = nmforderconsensus(consensus,kstart,kend);

clear ordcons clustid ordindex ; % Conserve memory.

cophCoefVector(1:kend)               = coph; % Values in coph(1:kstart-1) should be ignored
alreadyTried(kstart:kend)            = 1;    % No need to re-compute NMF at these kHat values.

% Parameters for NTU NMF code.
tol       = 1e-3;
timelimit = 5*24*60*60; % Allow up to five days (?!)
maxiter   = 128000;

% Initialize bracketing points 'lowVal' and 'hiVal'.  Determine starting stepSize.
loVal      = 2;
hiVal      = maxNC;
hiVal      = maxNC/2; % Don't try any values higher than half the rank.
%hiVal      = ceil(sqrt(maxNC)); % Don't try any values higher than sqrt(rank).
stepSize   = (hiVal-loVal)/subsamplingRate;

% Loop over possible values for number of components.
% Start indexing nc at '2' instead of '1'; otherwise nmf_taiwan code misbehaves.
finished = 0;

while ( finished == 0 )
    % Determine the indices for which FY and BICs will be computed this round.
    % Subsample the interval [lowVal, hiVal] by a factor 'stepSize'.
%subsamplingRate
%stepSize
%loVal
%keyboard
%if ( isempty(stepSize) | isempty(subsamplingRate) | isempty(loVal) )
%disp('EMPTY MATRIX FOUND!')
%    save junk.mat
%end
    indexList = round([0:subsamplingRate]*stepSize)+loVal;

    % Loop over the indices just computed, compute their Volumes.
    for nc=indexList
system(sprintf('touch /data/jmm97/NormalizationStudy/brunet.nc.%d.chk',nc));
if ( VERBOSE )
    disp(sprintf('mitEstimateK_ntuNMF_Multires: nc = %d\n',nc))
end
        if ( nc > length(alreadyTried) )
            continue
        end

        % If this value for 'k' hasn't been tried yet, evaluate it.
        if ( alreadyTried(nc) == 0 )

            % Brunet et al.'s cophenetic coefficient method.
            % Using the National Taiwan University NMF code to compute the NMF.
            kstart                               = nc;
            kend                                 = nc;
            consensus                            = nmfconsensus_taiwan(V',kstart,kend,numRuns,0);
            [ ordcons, clustid, ordindex, coph ] = nmforderconsensus(consensus,kstart,kend);
            clear ordcons clustid ordindex ; % Conserve memory.
            cophCoefVector(kend:kend)            = coph(kend:kend);

            % Update alreadyTried(nc).
            alreadyTried(nc) = 1;
        end; % if alreadyTried
    end; % for nc

    % Construct vector of all Volumes computed thus far, determine the Volume inflection point.
    alreadyTriedIndices = find(alreadyTried==1);
    cophCoefTmpList     = cophCoefVector(alreadyTriedIndices);

    % NEW: compute the cophenetic threshold from the maximum cophenetic
    % correlation computed thus far.
    maxCoph = max(cophCoefTmpList);
    cophThreshold = cophThreshFactor * maxCoph;

    % Update kHat, the current estimate for k.
    % Find the maximum index such that the cophenetic coefficient is at least
    % cophThreshold.
    indexOfIndex = max(find(cophCoefTmpList>=cophThreshold));
    kHat         = alreadyTriedIndices(indexOfIndex);

    % If we have reached the finest stepsize, we are finished.
    if ( stepSize <= 1 )
        finished = 1; 

    % Else, proceed to the next finest resolution in the WHILE loop by reducing the stepSize
    % and narrowing the search range [loVal,hiVal].
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

% Return current estimate of 'k', kHat.
return
