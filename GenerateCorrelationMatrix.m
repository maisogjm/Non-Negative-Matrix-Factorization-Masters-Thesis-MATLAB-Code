function [ C ] = GenerateCorrelationMatrix(nRows,offDiagVal)
% nRows      : number of rows (and number of columns, for a square matrix)
% offDiagVal : value for off-diagonal elements.

% Initialize C to identity matrix.
C = eye(nRows);

% Set off-diagonal elements.
C(find(C~=1)) = offDiagVal;

