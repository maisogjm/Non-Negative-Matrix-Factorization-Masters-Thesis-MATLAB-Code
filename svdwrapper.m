function [ u, s, v ] = svdwrapper( x, varargin )
%
% FORMAT [ returnVal ] = svdwrapper( x )
%
% FORMAT [ returnVal ] = svdwrapper( x, nu )
%
% FORMAT [ returnVal ] = svdwrapper( x, nu, nv )
%
% FORMAT [ returnVal ] = svdwrapper( x, nu, nv, verbose )
%
% The optional 4th argument 'verbose' is a boolean, set either to 'true' or 'false'.
%
% svd can crash.  Seems more common when cols outnumber rows
%
% Note: Art suggested trying the MATLAB function SVDS to compute reduced rank SVDs.
% So, this function is actually not called in my MATLAB translation.
% I have translated it to MATLAB anyway.

% Grab values out of 'vargin' for variables 'nu' and 'nv'.
% I suppose this section could use more error checking, like I did with 'verbose' below.
                                            % Below are line numbers and statements in the original R code.
if( nargin < 2 )                            % 6: if( missing(nu) )
    nu = min(size(x));                      % 6: nu = min(dim(x))
else
    nu = varargin{1};                       % 'nu' was given in the function call.
end
if( nargin < 3 )                            % 7: if( missing(nv) )
    nv = min(size(x))                       % 7: nv = min(dim(x))
else
    nv = varargin{2};                       % 'nv' was given in the function call.
end

% Set variable 'verbose' if an optional 4th argument was given.
if ( nargin == 4 )
    % Grab optional 4th argument, stuff it into variable 'verbose'.
    verbose = varargin{3};

    % Make sure 'verbose' is a logical value, and not something else like a struct.
    if ( class(verbose) ~= 'logical' )
        error('svdwrapper: ''verbose'' must be a logical scalar!')
    end
else
    verbose = true;                         % 4: verbose=T
end

% gotit = false; Not necessary w TRY/CATCH  % 9: gotit = F

try                                         % 11: try
    [ u, s, v ] = svd(x);                   % 11: svdx = svd(x,nu,nv);
    return

catch                                       % SVD failed, try svd of the transverse.
    try                                     % 15: try
        [ u, s, v ] = svd(x');              % 15: svdtx = svd(t(x),nv,nu);
        if ( verbose )                      % 26: if( verbose )
            % 26: print("svd(x) failed but svd(t(x)) worked.")
            disp('svd(x) failed but svd(t(x)) worked.');
        end
        temp = u;                           % 28: temp    = svdtx$u
        u    = v;                           % 29: svdtx$u = svdtx$v
        v    = temp;                        % 30: svdtx$v = temp
    catch                                   % 17: if( !gotit )
        nnan = sum(isnan(x(:)));            % 18: nnan = sum( is.nan(x) )
        if( nnan > 0 )                      % 19: if( nnan > 0 )
            error('x contained NaNs')       % 20: stop("x contained NaNs")
        else                                % 21: else{
            waybadx = x;                    % 22: waybadx <<- x
            % 23: stop("svd(x) and svd(t(x)) both failed without NaNs.\nExecuted waybadx <<- x.")
            error('svd(x) and svd(t(x)) both failed without NaNs. Executed waybadx = x.')
        end
    end; % Inner TRY
end; % Outer TRY

% As per MATLAB syntax, 'u', 'v', and 's' will be returned, now that they are set.
% svdtx
