function [ X ] = SubtractAbsoluteMinimum( X )
% This a very simple method to make positive data.
% The minimum negative value is substracted to every
% single cell of the data matrix.

X = X - min(X(:));
