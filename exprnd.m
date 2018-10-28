function [ expDev ] = exprnd(N, M, bbeta, varargin )
%
% FORMAT expDev = exprnd(N, M, bbeta)
%
% FORMAT expDev = exprnd(N, M, bbeta, ranSeed)
%
% N is the number of rows desired.
% M is the number of columns desired.
% bbeta is the exponential distribution parameter
% ranSeed is a random seed given as an optional 4th input argument
%
%
% The output expDev is an NxM matrix of pseudorandom numbers drawn from an exponential distribution.

% Make sure N and M are both greater than zero.
if ( N <= 0 )
    error('exprnd: ''N'' must be a positive number!')
end
if ( M <= 0 )
    error('exprnd: ''M'' must be a positive number!')
end

% Make sure bbeta is greater than zero.
if ( bbeta <= 0 )
    error('exprnd: ''bbeta'' must be a positive number!')
end

% Set random seed if an optional 4th argument was given.
if ( nargin == 4 )
    % Grab optional 4th argument, set random seed.
    ranSeed = varargin{1};

    % Make sure 'check' is a double-precision value, and not something else like a struct.
    if ( class(ranSeed) == 'double' )
        rand('state',ranSeed); % Set the random number seed.
    else
        error('exprnd: ''ranSeed'' must be a double-precision scalar!')
    end
end

% Generate a vector of uniformly distributed random numbers, none of which is zero.
expDev      = rand(N,M);
zeroIndices = find(expDev == 0 );
while( ~isempty(zeroIndices) )
    N                   = length(zeroIndices); % If there are any zero values...
    expDev(zeroIndices) = rand(N,1);           % fill them in with random numbers...
    zeroIndices         = find(expDev == 0 );  % until there are no zero values left.
end

% Transform uniform deviates to exponential deviates.
% Source: Press WH, Teukolsky SA, Vetterling WT, Flannery BP,
% Numerical Recipes, 3rd edition, Cambridge:Cambridge University Press (2009), p. 362.
expDev = -log(expDev)/bbeta;
