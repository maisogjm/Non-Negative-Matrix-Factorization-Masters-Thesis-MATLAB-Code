function [ X ] = FoldDataByRows( X )
% This approach was used by Kim and Tidor (Genome Res. 2003) for the
% analysis of log-transformed gene expression data. Every row (item)
% is represented in two new rows of a new matrix. The first one is
% used to indicate positive expression (up-regulation) and the second
% one to indicate a negative expression value (down-regulation).
% This process doubles the number of rows of the data set.

% If there are no negative values in X, return X.
if ( sum(X(:)<0) == 0 )
    return;
end

% Append a negative version of X to X, doubling the number of rows.
X = [ X ; -X ];

% Set negative values to 0.
X(find(X<0)) = 0;
