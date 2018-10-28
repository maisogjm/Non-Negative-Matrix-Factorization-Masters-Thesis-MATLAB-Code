function W = StepFunctionColumns(nRows, nCols)
%
% This function generates a matrix whose rows are simple orthogonal step functions
%
% INPUTS:
% nRows  - desired number of rows
% nCols  - desired number of columns
% 
% OUTPUT:
%
% W - matrix whose rows are simple orthogonal step functions

% If nRows < nCols, exit out with error.
if ( nRows < nCols )
    error('StepFunctionRows: nRows must be greater than or equal to nCols!')
end

W = zeros(nRows,nCols);

% Loop over columns, design step functions
stepFactor = nRows / nCols ;
startIndexTmp = 1; % Initialize startIndex.
for i=1:nCols
    % Compute endIndex from startIndex.
    endIndexTmp = startIndexTmp + stepFactor - 1;
    startIndex  = round(startIndexTmp);
    endIndex    = round(endIndexTmp);

    % If this is the last column, make sure the 1's should go all the way to the last row.
    if ( i == nCols)
        endIndex = nRows;
    end

    % Set select elements in the i-th column of W to 1.
    W(startIndex:endIndex,i) = ones(endIndex-startIndex+1,1);

    % Update startIndex.
    startIndexTmp = endIndexTmp + 1;
end
