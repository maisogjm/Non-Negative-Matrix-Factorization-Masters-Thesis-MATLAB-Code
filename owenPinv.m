% This function was originally named 'pinv', but there is already a MATLAB function named PINV.
% To avoid "namespace collisions", I have renamed this function to OWENPINV.  By MATLAB rules,
% the filename must also be changed accordingly.
function [ xInv ] = owenPinv( x, varargin )
%
% FORMAT function [ xInv ] = owenPinv(x)
%
% FORMAT function [ xInv ] = owenPinv(x, k)
%
% FORMAT function [ xInv ] = owenPinv(x, k, eta)
%
% FORMAT function [ xInv ] = owenPinv(x, k, eta, eps)
%
% function [ xInv ] = owenPinv( x,k,eta=1e-100,eps=1e-9)
% Pseudo inverse of x, optionally truncated to rank k

% Obtain variables from 'varargin'.           % Below are line numbers and statements in the original R code.
if ( nargin >= 2 )                            % 5: if( missing(k)    )
    k = varargin{1};
else
    k = min(size(x));                         % 5: k = min(dim(x))
end

if ( nargin >= 3 )
    eta = varargin{2};
else
    eta = 1e-100;                             % 3: eta=1e-100
end

% The original R code had eps defaulting to a value of 1e-9.
% In this MATLAB translation, I'll let it default to the pre-defined MATLAB value for eps.
% In MATLAB Version 7.1.0.183 (R14), this value is 2.2204e-16.
if ( nargin >= 4 )
    epsLocal = varargin{3};
else
    epsLocal = eps;                           % 3: eps=1e-9
end

if ( k > min(size(x)) )                       % 6: if( k>min(dim(x)) )
    k = min(size(x));                         % 6: k = min(dim(x))
end

% Art suggested trying the MATLAB function SVDS to compute reduced rank SVDs.
[ u, d, v ] = svds(x,k);                      % 8: svdx = svdwrapper(x,k,k)
d           = diag(d);                        % MATLAB returns a diagonal matrix, whereas R returns a vector.
                                              % Force 'd' to be a vector just to keep the code similar to the original.

thresh = epsLocal*max(d) + eta;               % 10: thresh = eps*max(svdx$d) + eta

d = d(1:k);                                   % 12: svdx$d = svdx$d[1:k]

dinv = d;                                     % 14: dinv = svdx$d
dinv = ones(length(d),1) ./ dinv;             % 15: dinv = 1.0/dinv
dinv( find( abs(dinv) >= 1/thresh ) ) = 0;    % 16: dinv[abs(dinv)>=1/thresh] = 0 % dinv  < 0 is possible numerically

for j=1:k;                                    % 18: for( j in 1:k )
    v(:,j) = v(:,j) * dinv(j);                % 18: svdx$v[,j] = svdx$v[,j]*dinv[j]
end

% Compute pseudo-inverse.
xInv = v * u';                                % 20: svdx$v %*% t(svdx$u)

% As per MATLAB syntax, xInv will be returned, now that it is set.

%=======================================================================
%=======================================================================
% The LRA function has originally in the file pinv.R
% However, only the function GABR invokes LRA.
% So, it seemed best to make LRA a local subfunction in the file gabr.m.
%
% ## rank k approximation to x
% ## we don't let k be missing
%
% lra = function( x,k ){
%
%  if( k > min(dim(x)) )
%    return(x)
%
%  svdx = svdwrapper(x,k,k)
%  d    = svdx$d
%  for( j in 1:min(k,length(d)) )svdx$v[,j] = svdx$v[,j]*d[j]
%
%  svdx$u[,1:k] %*% t(svdx$v[,1:k])
% }

%=======================================================================
%=======================================================================
% R has a function RANGE that is invokved by PINVTESTS.
% We'll have to define one for MATLAB.
function [ minX, maxX ] = range(X)
    minX = min(X(:));
    maxX = max(X(:));
return

%=======================================================================
%=======================================================================
function [ returnVal ] = pinvtests(varargin)
%
% FORMAT [ returnVal ] = pinvtests() 
%
% FORMAT [ returnVal ] = pinvtests(n) 
%
% FORMAT [ returnVal ] = pinvtests(n, k) 
%
% FORMAT [ returnVal ] = pinvtests(n, k, rep) 
%
% FORMAT [ returnVal ] = pinvtests(n, k, rep, seed) 
%
% Test four pseudo-inverse conditions on random matrices
%
%function [ returnVal ] = pinvtests(n = 10, k=5, rep=10, seed )

if ( nargin >= 4 )                            % 44: if( !missing(seed) )
    seed = varargin{4};
    rand('state',seed);                       % 44: set.seed(seed)
end

for r=1:rep                                   % for( r in 1:rep )
    x = rand(n,k);                            % 51: x    = matrix( runif(n*k),n,k )
    px = owenPinv(x);                         % 52: px   = pinv(x)

    i    = range(  x * px * x  -  x );        % 54: i    = range(  x %*% px %*% x  -  x )
    ii   = range( px *  x * px - px );        % 55: i   = range( px %*%  x %*% px - px ) 
    iii  = range( (x*px)' - x*px );           % 56: iii  = range( t(x%*%px) - x%*%px )
    iv   = range( (px*x)' - px*x );           % 57: iv   = range( t(px*x) - px*x );

    disp([ i ii iii iv ]);                    % 59: print(c(i,ii,iii,iv))
end
%
% Test on matrix of zeros too
%
x = zeros(n,k);                               % 65: x = matrix(0,n,k)
disp(range(owenPinv(x)));                     % 66: print(range(pinv(x)))

return
