function [ rowgps, colgps ] = genbicvgps(m, n, k, l )
%
% FORMAT [ rowgps, colgps ] = genbicvgps(m, n, k, l)
%
% FORMAT [ rowgps, colgps ] = genbicvgps(m, n, k, l, seed)
%
% Generate groups for (k,l) bi-cross-validation
%   Split rows  1:m into k groups of nearly equal size
%   and columns 1:n into l groups of nearly equal size
%
% If an optional 5th scalar argument 'seed' is given, it will be used to re-seed the
% random number generator.
%
% Code is not optimized for large k or l
%
% Row groups with small index appear ceiling(m/k) times
% Others (if any) appear ceiling(m/k)-1 times

% Round 'k' and 'l' to nearest integer.
                                   % Below are line numbers and statements in the original R code.
k = round(k) ; l = round(l) ;      % 14: k = round(k) ; l = round(l)

% Set RNG seed if an optional 5th argument was given.
% if ( nargin == 2 )                 % 16: if( !missing(seed) )
    % Grab optional 5th argument, stuff it into variable 'seed'.
%     seed = varargin{1};

    % Make sure 'seed' is a double-precision value, and not something else like a struct.
%     if ( class(seed) == 'double' )
%         rand('state',seed);        % 16: set.seed(seed)
%     else
%         error('genbicvgps: ''seed'' must be a double precision scalar!')
%     end
% end

% Rows.  The MATLAB equivalent of the R function SAMPLE is the function RANDPERM.
rowBase   = 1:k;                    % 18: 1:k
N      = ceil(m/k);
rowgps = zeros(1,k*N);              % Pre-allocate memory.
for p = 1:N                         % 18: ceiling(m/k)
    start  = (p-1)*k+1;
    finish = (p-1)*k+k;
    rowgps(start:finish) = rowBase; % 18: rowgps = rep( 1:k, ceiling(m/k) )
end
rowgps = rowgps(1:m);               % 19: rowgps = rowgps[1:m]
rowgps = rowgps(randperm(m));       % 20: rowgps = rowgps[sample(m)]

% Columns.  The MATLAB equivalent of the R function SAMPLE is the function RANDPERM.
colBase   = 1:l;                    % 22: 1:l
M      = ceil(n/l);
colgps = zeros(1,l*M);              % Pre-allocate memory.
for p = 1:M                         % 22: ceiling(n/l)
    start  = (p-1)*l+1;
    finish = (p-1)*l+l;
    colgps(start:finish) = colBase; % 18: rowgps = rep( 1:k, ceiling(m/k) )
end
colgps = colgps(1:n);               % 23: colgps = colgps[1:n]
colgps = colgps(randperm(n));       % 24: colgps = colgps[sample(n)]

% As per MATLAB syntax, rowgps and colgps will be returned, now that they are set.
% list( rowgps=rowgps, colgps=colgps )
