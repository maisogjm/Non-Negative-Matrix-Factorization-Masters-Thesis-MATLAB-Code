function [ X ] = ScaleColumnsThenNormalizeRows( X )
% This is the approach proposed by Getz, et al. (PNAS 2000) that
% first divide each column by its mean and then normalize each row.
% Using suggested method to avoid FOR loops found here:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/155998

% First divide each column by its mean.
[ mRows nCols ] = size(X);  % Get dimensions of X.
M = mean(X,1);              % Compute mean of each column.
M = M(ones(mRows,1),:);     % Form matrix with duplicate values within column.
X = X./M;                   % Scale each column by its mean.

% Then normalize each row ("such that its mean vanishes and its norm is one").
X  = SubtractMeanByRows(X); % Set row means to 0.
SS = sqrt(sum(X.^2,2));     % Compute row norms.
SS = SS(:,ones(1,nCols));   % Form matrix with duplicate values within row.
X  = X./SS;                 % Set row norm to 1.
