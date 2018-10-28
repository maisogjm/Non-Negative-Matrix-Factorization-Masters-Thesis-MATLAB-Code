function [ X ] = ExponentialScale( X )
% Data is exponentially scaled to make it positive.
% This is an inverse operation of a logarithmic transformation.

X = exp(X);
