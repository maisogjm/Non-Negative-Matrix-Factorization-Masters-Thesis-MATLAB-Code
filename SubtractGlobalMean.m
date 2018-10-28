function [ X ] = SubtractGlobalMean( X )
% The global mean of the data matrix is calculated
% and then subtracted from all data items..

X = X - mean(X(:)); % Subtract global mean.
