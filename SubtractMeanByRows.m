function [ X ] = SubtractMeanByRows( X )
% The mean for each row of the data matrix is calculated
% and then substracted from all data items of that row.
% Using suggested method to avoid FOR loops found here:
% http://www.mathworks.com/matlabcentral/newsreader/view_thread/155998

[ mRows nCols ] = size(X); % Get dimensions of X.
V = mean(X,2);             % Compute mean of each row.
V = V(:,ones(1,nCols));    % Form matrix with duplicate values within row.
X = X-V;                   % Normalize X.
