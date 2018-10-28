function [ simdat, C ] = SimulateData_KD16()
%
% FORMAT [ simdat, C ] = SimulateData_KD16()
%
% Simulate data with 16 clusters, using Karthik's style.
% 50-block genes will be used to define clusters.
% Also returns cross-correlation matrix of clusters.

global numClusters
numClusters = 16;

% Initialize 1000x60 matrix with noise
% ~ EXp(Beta=1) distribution
simdat = exprnd(1000,60,1);

% Initialize matriX for computing the correlation.
X = zeros(1000,numClusters);

% Construct cluster #1.
[ simdat, X ] = SetBlock(simdat,X,1,1,1/40);
[ simdat, X ] = SetBlock(simdat,X,1,2,1/1);
[ simdat, X ] = SetBlock(simdat,X,1,3,1/1);
[ simdat, X ] = SetBlock(simdat,X,1,4,1/1);

% Construct cluster #2.
[ simdat, X ] = SetBlock(simdat,X,2,1,1/80);
[ simdat, X ] = SetBlock(simdat,X,2,2,1/80);
[ simdat, X ] = SetBlock(simdat,X,2,3,1/1);
[ simdat, X ] = SetBlock(simdat,X,2,4,1/1);

% Construct cluster #3.
[ simdat, X ] = SetBlock(simdat,X,3,1,1/1);
[ simdat, X ] = SetBlock(simdat,X,3,2,1/120);
[ simdat, X ] = SetBlock(simdat,X,3,3,1/1);
[ simdat, X ] = SetBlock(simdat,X,3,4,1/1);

% Construct cluster #4.
[ simdat, X ] = SetBlock(simdat,X,4,1,1/40);
[ simdat, X ] = SetBlock(simdat,X,4,2,1/80);
[ simdat, X ] = SetBlock(simdat,X,4,3,1/1);
[ simdat, X ] = SetBlock(simdat,X,4,4,1/120);

% Construct cluster #5.
[ simdat, X ] = SetBlock(simdat,X,5,1,1/80);
[ simdat, X ] = SetBlock(simdat,X,5,2,1/80);
[ simdat, X ] = SetBlock(simdat,X,5,3,1/40);
[ simdat, X ] = SetBlock(simdat,X,5,4,1/120);

% Construct cluster #6.
[ simdat, X ] = SetBlock(simdat,X,6,1,1/80);
[ simdat, X ] = SetBlock(simdat,X,6,2,1/80);
[ simdat, X ] = SetBlock(simdat,X,6,3,1/120);
[ simdat, X ] = SetBlock(simdat,X,6,4,1/120);

% Construct cluster #7.
[ simdat, X ] = SetBlock(simdat,X,7,1,1/40);
[ simdat, X ] = SetBlock(simdat,X,7,2,1/120);
[ simdat, X ] = SetBlock(simdat,X,7,3,1/120);
[ simdat, X ] = SetBlock(simdat,X,7,4,1/40);

% Construct cluster #8.
[ simdat, X ] = SetBlock(simdat,X,8,1,1/1);
[ simdat, X ] = SetBlock(simdat,X,8,2,1/1);
[ simdat, X ] = SetBlock(simdat,X,8,3,1/80);
[ simdat, X ] = SetBlock(simdat,X,8,4,1/120);

% Construct cluster #9.
[ simdat, X ] = SetBlock(simdat,X,9,1,1/120);
[ simdat, X ] = SetBlock(simdat,X,9,2,1/40);
[ simdat, X ] = SetBlock(simdat,X,9,3,1/120);
[ simdat, X ] = SetBlock(simdat,X,9,4,1/40);

% Construct cluster #10.
[ simdat, X ] = SetBlock(simdat,X,10,1,1/40);
[ simdat, X ] = SetBlock(simdat,X,10,2,1/80);
[ simdat, X ] = SetBlock(simdat,X,10,3,1/1);
[ simdat, X ] = SetBlock(simdat,X,10,4,1/120);

% Construct cluster #11.
[ simdat, X ] = SetBlock(simdat,X,11,1,1/120);
[ simdat, X ] = SetBlock(simdat,X,11,2,1/1);
[ simdat, X ] = SetBlock(simdat,X,11,3,1/120);
[ simdat, X ] = SetBlock(simdat,X,11,4,1/120);

% Construct cluster #12.
[ simdat, X ] = SetBlock(simdat,X,12,1,1/1);
[ simdat, X ] = SetBlock(simdat,X,12,2,1/120);
[ simdat, X ] = SetBlock(simdat,X,12,3,1/80);
[ simdat, X ] = SetBlock(simdat,X,12,4,1/80);

% Construct cluster #13.
[ simdat, X ] = SetBlock(simdat,X,13,1,1/120);
[ simdat, X ] = SetBlock(simdat,X,13,2,1/120);
[ simdat, X ] = SetBlock(simdat,X,13,3,1/120);
[ simdat, X ] = SetBlock(simdat,X,13,4,1/1);

% Construct cluster #14.
[ simdat, X ] = SetBlock(simdat,X,14,1,1/1);
[ simdat, X ] = SetBlock(simdat,X,14,2,1/120);
[ simdat, X ] = SetBlock(simdat,X,14,3,1/40);
[ simdat, X ] = SetBlock(simdat,X,14,4,1/120);

% Construct cluster #15.
[ simdat, X ] = SetBlock(simdat,X,15,1,1/80);
[ simdat, X ] = SetBlock(simdat,X,15,2,1/1);
[ simdat, X ] = SetBlock(simdat,X,15,3,1/120);
[ simdat, X ] = SetBlock(simdat,X,15,4,1/120);

% Construct cluster #16.
[ simdat, X ] = SetBlock(simdat,X,16,1,1/40);
[ simdat, X ] = SetBlock(simdat,X,16,2,1/1);
[ simdat, X ] = SetBlock(simdat,X,16,3,1/1);
[ simdat, X ] = SetBlock(simdat,X,16,4,1/80);

C = corr(X);

%--------------------------------------------------
% Subfunctions for setting matrix element values
% by submatrix (a block of genes).
%--------------------------------------------------
function [ M, X ] = SetBlock(M,X,cluster,block,val)
    % Some constants.
    global numClusters
    blockRows   = 50; % Number of rows per block.
    [ nRows, nCols ] = size(M); % Get dimensions;

    % Define submatrix boundaries.
    rowStart  = ( blockRows * (block-1) ) + 1;
    rowEnd    = rowStart + blockRows - 1;
    colFactor = nCols / numClusters; % Length factor.
    colStart  = round((cluster-1)*colFactor)+1;
    colEnd    = round(cluster*colFactor);
    blockCols = colEnd - colStart + 1;

    % Set submatriX.
    M(rowStart:rowEnd,colStart:colEnd) = ...
        exprnd(blockRows,blockCols,val);
    
    % Set vector.  This will be used to compute
    % the cross-correlation matriX.
    X(rowStart:rowEnd,cluster) = (1/val)*ones(blockRows,1);
return
