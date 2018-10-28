% Make random data.
m        = 100;  % Number of rows.
n        = 1000; % Number of columns.
numComps = 2;    % Number of components.
%x        = ParkKimCichokiHybridSimData(m,n,numComps,0.4,0.05);
x        = rand(m,n);

% Generate groups for (k,l) bi-cross-validation
k = 3;
l = 30;
[ rowgps, colgps ] = genbicvgps(m, n, k, l);

% Take those groups and do the BCV using function gabr
ks =  1:4;
[ gI ] = gengabcvsvd( x,  [ rowgps, colgps ], ks )
