function result = RSSQ(X,Xhat)
% Function returns a norm of the difference matrix: the relative sum of square differences.
% This norm is computed as the sum of the squared differences,
% divided by the sum of the squared values of the original matrix.
diffMat = X - Xhat;
numer   = sum(sum(diffMat.*diffMat));
denom   = sum(sum(X.*X));
if ( denom ~= 0 )
    result = numer / denom;
else
    disp('rssq: denominator is zero, which implies X is all zero!')
    result = NaN;
end
