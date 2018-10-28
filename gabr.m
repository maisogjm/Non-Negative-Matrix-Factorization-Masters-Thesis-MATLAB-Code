function [ A ] = gabr(x, k, iout, jout, varargin)
%
% FORMAT gabr = function(x, k, iout, jout)
%
% FORMAT gabr = function(x, k, iout, jout, check)
%
% FORMAT gabr = function(x, k, iout, jout, check, versionNumber)
%
% Gabriel residual, versions I and II
%
% X is an m by n matrix
% iout is a list of r distinct values from 1:m
% jout is a list of s distinct values from 1:n
%
% gabr = function( x, k, iout, jout, check=T, version="I" )

% Set variable 'check' if an optional 5th argument was given.
if ( nargin == 5 )
    % Grab optional 5th argument, stuff it into variable 'verbose'.
    verbose = varargin{1};

    % Make sure 'check' is a logical value, and not something else like a struct.
    if ( class(check) ~= 'logical' )
        error('gabr: ''check'' must be a logical scalar!')
    end
else;                                         % Below are line numbers and statements in the original R code.
    check = true;                             % 7: check=T
end

% Set variable 'versionNumber' if an optional 6th argument was given.
if ( nargin == 6 )
    % Grab optional 5th argument, stuff it into variable 'versionNumber'.
    versionNumber = varargin{2};

    % Make sure 'versionNumber' is a logical value, and not something else like a struct.
    if ( class(versionNumber) ~= 'char' )
        error('gabr: ''versionNumber'' must be a character string!')
        if ( ( versionNumber ~= 'I' ) && ( versionNumber ~= 'II' ) )
            error('gabr: ''versionNumber'' must be either "I" or "II"!')
        end
    end
else
    versionNumber = 'I';                      % 7: version="I"
end

if ( check )                                     % 9: if( check ){
    if ( ~isnumeric(x) )                         % 10: if( !is.matrix(x)       )
        error('x is not numeric')
    else
        if ( isempty(x) )
            error('x is empty')
        else
            [ nRows, nCols ] = size(x);
            if ( ( nRows == 1 ) || ( nCols == 1 ) )
                error('x is not a matrix');      % 10: error('x is not a matrix')
            end
        end
    end
    if ( max(iout(:)) > nRows )                  % 11: if( max(iout) > nrow(x) )
        error('max iout too large')              % 11: stop("max iout too large")
    end
    if ( min(iout(:)) < 1 )                      % 12: if( min(iout) < 1       )
        error('min iout too small')              % 12: stop("max iout too large")
    end
    if( max(jout(:)) > nCols )                   % 13: if( max(jout) > ncol(x) )
        error('max jout too large')              % 13: stop("max jout too large")
    end
    if ( min(jout(:)) < 1 )                      % 14: if( min(jout) < 1       )
        error('min jout too small')              % 14: stop("min jout too small")
    end

    if( length(unique(iout)) < length(iout) )    % 16: if(  length( unique(iout)) < length(iout) )
        error('duplicates in iout')              % 16: stop("duplicates in iout")
    end
    if( length(unique(jout)) < length(jout) )    % 17: if(  length( unique(jout)) < length(jout) )
        error('duplicates in jout')              % 17: stop("duplicates in jout")
    end

    if( k<0 )                                    % 19: if( k<0 )
        error('negative k')                      % 19: stop("negative k")
    end
end; % IF check                                  % 20: }
if ( k > max(size(x)) )                          % 21: if( k>max(dim(x)))
    k=max(size(x))                               % 21: k=max(dim(x))
end

A            = x(iout,jout);                     % 23: A = x[  iout,  jout, drop=F ]

xTmp         = x;
xTmp(:,jout) = [];                               % Drop column jout.
B            = xTmp(iout, :);                    % 24: B = x[  iout, -jout, drop=F ]

xTmp         = x;
xTmp(iout,:) = [];                               % Drop row iout.
C            = xTmp(:,jout);                     % 25: C = x[ -iout,  jout, drop=F ]

xTmp         = x;
xTmp(iout,:) = [];
xTmp(:,jout) = [];                               % Drop column jout and row iout.
D            = xTmp;                             % 26: D = x[ -iout, -jout, drop=F ]

if ( k==0 )                                      % 28: if( k==0 )
    return;                                      % 28: return(A)
end

% print("got abcd")
% print(A);  print(B);  print(C);  print(D)

if ( versionNumber == 'II' )                           % 33: if( version == "II" ){
    B = lra(B,k);                                % 34: B = lra(B,k)
    C = lra(C,k);                                % 35: C = lra(C,k)
end                                              % 36: }

% A = A - B * owenPinv(D,k) * C;                   % 38: A - B %*% pinv(D,k) %*% C
A = A - B * pinv(D,k) * C;                   % 38: A - B %*% pinv(D,k) %*% C

return

%=======================================================================
%=======================================================================
function [] = gabrtest( varargin )
%
% FORMAT
%
% gabrtest = function( m = 10, n=5, r=2, s=2, k=3, rep=10, seed )

% Obtain variables from 'varargin'.
if ( nargin >=1 )
    m = varargin{1};
else
    m = 10;                                      % 42: m = 10
end

if ( nargin >= 2 )
    n = varargin{2};
else
    n = 5;                                       % 42: n=5
end

if ( nargin >= 3 )
    r = varargin{3};
else
    r = 2;                                       % 42: r=2
end

if ( nargin >= 4 )
    s = varargin{4};
else
    s = 2;                                       % 42: s=2
end

if ( nargin >= 5 )
    k = varargin{5};
else
    k = 3;                                       % 42: k=3
end

if ( nargin >= 6 )
    rep = varargin{6};
else
    rep = 10;                                    % 42: rep=10
end

if ( nargin >= 7 )                               % 44: if( !missing(seed) )
    seed = varargin{7};
    rand('state',seed);                          % 44: set.seed(seed)
end

% spike matrix
x      = zeros(m,n);                             % 47: x = matrix(0,m,n)
x(1,1) = 1;                                      % 48: x[1,1] = 1

iout = 1:r;                                      % 50: iout = 1:r
jout = 1:s;                                      % 51: jout = 1:s

disp(gabr(x,k,iout,jout,T,'I' ));                % 53: print( cbind(gabr(x,k,iout,jout,T,"I"),
disp(gabr(x,k,iout,jout,T,'II'));                % 54: gabr(x,k,iout,jout,T,"II")))

% stripe matrix
x      = zeros(m,n);                             % 57: x = matrix(0,m,n)
x(1,:) = 1;                                      % 58: x[1,] = 1

iout = 1:r;                                      % 60: iout = 1:r
jout = 1:s;                                      % 61: jout = 1:s

% angle matrix
x      = zeros(m,n);                             % 64: x = matrix(0,m,n)
x(1,:) = 1;                                      % 65: x[1,] = 1
x(:,1) = 1;                                      % 66: x[,1] = 1

iout = 1:r;                                      % 68: iout = 1:r
jout = 1:s;                                      % 69: jout = 1:s

disp(gabr(x,k,iout,jout,T,'I' ));                % 71: print( cbind(gabr(x,k,iout,jout,T,"I"),
disp(gabr(x,k,iout,jout,T,'II'));                % 72: gabr(x,k,iout,jout,T,"II")))

for i = 1:rep                                    % 74: for( i in 1:rep ){
    x = round( matrix(runif(m*n),m,n), 3 )       % 75: x = round( matrix(runif(m*n),m,n), 3 )
    % 76: print( gabr(x,k,iout,jout,T,"I")-gabr(x,k,iout,jout,T,"II"))
    disp( gabr(x,k,iout,jout,T,'I')-gabr(x,k,iout,jout,T,'II'))
end; % FOR i                                     % 77: }
return

%=======================================================================
%=======================================================================
% The following function LRA was originally in the file pinv.R.
% However, only the function GABR invokes LRA.
% So, it seemed best to make LRA a local subfunction in the file gabr.m.
function [ returnVal ] = lra( x,k )
%
% FORMAT [ returnVal ] = lra( x,k )
%
% rank k approximation to x
% we don't let k be missing

if( k > min(size(x)) )                           % 29: if( k > min(dim(x)) )
    returnVal = x;                               % 30: return(x)
    return;                                      % 30: return(x)
end

% Art suggested trying the MATLAB function SVDS to compute reduced rank SVDs..
[ u, d, v ] = svds(x,k);                         % 32: svdx = svdwrapper(x,k,k)
d           = diag(d);                           % 33: d    = svdx$d
for j = 1:min(k,length(d))                       % 34: for( j in 1:min(k,length(d)) )
    v(:,j) = v(:,j) * d(j);                      % 34: svdx$v[,j] = svdx$v[,j]*d[j]
end

returnVal = u(:,1:k) * v(:,1:k)';                % 36: svdx$u[,1:k] * t(svdx$v[,1:k])
return

