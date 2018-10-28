function result = RRSSQ(X,Xhat)
% Function returns a norm of the difference matrix: the relative root of the sum of squared differences.
% This norm is computed as the sum of the squared differences,
% divided by the sum of the squared values of the original matrix.
% Then the square root is taken of the result.
result = sqrt(RSSQ(X,Xhat));
