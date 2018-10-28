function Y = inconsistent(Z,depth)
%INCONSISTENT Inconsistent values of a cluster tree.
%   Y = INCONSISTENT(Z) computes the inconsistent values of the
%   cluster tree given by Z. Z is a M-1 by 3 matrix generated from the
%   function LINKAGE. 
%
%   Y = INCONSISTENT(Z, DEPTH) computes the inconsistent values
%   with depth equal to DEPTH.
%
%   The output Y is a M-1 by 4 matrix. Suppose S is the pruned
%   subtree starting from node i to a depth of DEPTH, then
%
%      Y(i,1) = Average distance between nodes in S
%      Y(i,2) = Standard deviation of distances between nodes in S
%      Y(i,3) = Number of nodes in S
%      Y(i,4) = (Distance of node i - Y(i,1))/Y(i,2)
%
%   When DEPTH is not given, the default value is set to be 2.
%
%   See also PDIST, LINKAGE, COPHENET, DENDROGRAM, CLUSTER, CLUSTERDATA.

%   Copyright 1993-2002 The MathWorks, Inc. 
%   $Revision: 1.10 $

if nargin < 2, depth = 2; end

m = size(Z,1);

Y = zeros(m,4);

for k = 1:m
   s = zeros(4,1);
   s = tracetree(Z, s, k, depth);
   Y(k,1) = s(1)/s(3); % average edge length 
   Y(k,2) = sqrt((s(2) - (s(1)*s(1))/s(3))/(s(3)-(s(3)~=1))); % standard deviation
   Y(k,3) = s(3); % number of edges 
   if Y(k,2) > 0
      Y(k,4) = (Z(k,3) - Y(k,1))/Y(k,2);
   else 
      Y(k,4) = 0;  
   end
end

function s = tracetree(Z,s,k,depth)
% iterative function to search down the tree.
m = size(Z,1)+1;
klist = zeros(m,1);
klist(1) = k;
dlist(1) = depth;
topk = 1;
currk = 1;
while(currk <= topk)
   k = klist(currk);
   depth = dlist(currk);
   s(1) = s(1) + Z(k,3); % sum of the edge lengths so far
   s(2) = s(2) + Z(k,3)*Z(k,3); % sum of the square of the edge length 
   s(3) = s(3) + 1; % number of the edges so far 

   if depth > 1 % depth is greater than 0, need to go down further 
      for i = Z(k,1:2);   % left and right subtree indices
         if i > m % node i is not a leaf, it has subtrees 
            topk = topk+1;
            klist(topk) = i-m;
            dlist(topk) = depth-1;
         end
      end
   end
   currk = currk+1;
end
