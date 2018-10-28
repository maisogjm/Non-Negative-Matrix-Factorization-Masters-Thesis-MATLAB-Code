function [ X ] = SubtractMeanByRowsAndThenByColumns( X )
% The mean for each row of the data matrix is calculated and then
% substracted from all data items of that row. In a subsequent
% step, the mean for each column of the data matrix is calculated
% and then substracted from all data items of that column.

X = SubtractMeanByRows(X);    % Normalize rows of X.
X = SubtractMeanByColumns(X); % Normalize columns of X.
