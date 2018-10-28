function [ gI ] = gengabcvsvd( x, rowgps, colgps, ks, verbose )
%
% FORMAT [ gI ] = gengabcvsvd( x,  rowgps, colgps, ks)
%
% FORMAT [ gI ] = gengabcvsvd( x,  rowgps, colgps, ks, verbose )
%
% Switch NNMF code to SVD
% may be extremely inefficient
% does not reuse svds
%
% [ gI ] = gengabcvsvd( x,  gps, ks, verbose=1 )

% Obtain variables from 'varargin'.      % Below are line numbers and statements in the original R code.
% if ( nargin >= 4 )
%     verbose = varargin{1};
% else
%     verbose = 1;                      % 5: verbose=1
% end

gI   = ks*0;                             % 7: gI   = ks*0

% rowgps = gps(1);                       % 9:  rowgps = gps$rowgps
% colgps = gps(2);                       % 10: colgps = gps$colgps

nr = max( rowgps );                      % 12: nr = max( rowgps )
nc = max( colgps );                      % 13: nc = max( colgps )

for l = 1:length(ks)                     % 15: for( l in 1:length(ks) ){
    for i = 1:nr                         % 16: for( i in 1:nr )
        for j = 1:nc                     % 16: for( j in 1:nc ){
            % The MATLAB command FIND is the equivalent to the R command WHICH.
            ro = find(rowgps==i);        % 17: ro = which(rowgps==i) % rows left out
            co = find(colgps==j);        % 18: co = which(colgps==j) % cols left out

            tmpMat = gabr(x, ks(l), ro, co, true, 'I' ); % 20: gabr(    x, ks[l], ro, co, T, "I" )

            % 20: gI[l]  = gI[l]  + sum(gabr(    x, ks[l], ro, co, T, "I" )^2)
            gI(l) = gI(l) + sum(sum(tmpMat .* tmpMat));

            if( verbose>0)               % 21: if( verbose>1)
                disp( [ ks(l), gI(l) ] ) % 21: print( c(ks[l],gI[l]) )
            end
        end; % FOR j

        if(verbose>0)                    % 24: if(verbose>0)
            disp( [ ks(l), gI(l) ] )     % 24: print( c(ks[l],gI[l]) )
        end
    end; % FOR i                         % 22: }
end; % FOR l                             % 25: }

