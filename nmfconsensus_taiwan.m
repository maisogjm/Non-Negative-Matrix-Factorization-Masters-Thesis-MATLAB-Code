function consensus = nmfconsensus_taiwan(a,kstart,kend,nloop,verbose)
%
% Jean-Philippe Brunet
% Cancer Genomics
% The Broad Institute
% brunet@broad.mit.edu
%
% This software and its documentation are copyright 2004 by the
% Broad Institute/Massachusetts Institute of Technology. All rights are reserved.
% This software is supplied without any warranty or guaranteed support whatsoever.
% Neither the Broad Institute nor MIT can not be responsible for its use, misuse,
% or functionality.
%
% Model selection for NMF
%
% a (n,m) : N (genes) x M (samples) original matrix
%           numerical data only. Must be positive.
%
% kstart, kend : range of values of k to test consensus for.
%
% nloop : number of initial conditions per k
%         (start with 10 or 20, check with more)
%
% verbose : prints iteration count and changes in connectivity matrix elements
%           if not set to 0
%
% consensus : 3d array of consensus matrices
%             dimensions : kend x M x M
%             Values in consensus(1:kstart-1,:,:) should be ignored
%

% test for negative values in v
if min(min(a)) < 0
    error('matrix entries can not be negative');
    return
end
if min(sum(a,2)) == 0
    error('not all entries in a row can be zero');
    return
end

[n,m]=size(a);

%disp(sprintf('nmfconsensus_taiwan: m = %d',m))
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Rows.%d.Columns.%d.chk',n,m))
%exit
consensus = zeros(kend,m,m);
conn      = zeros(m,m);
incr      = 0;

% Parameters for NTU code.
tol       = 1e-3;
timelimit = 5*24*60*60; % Allow up to five days (?!)
maxiter   = 128000;

for j=kstart:kend
    if verbose fprintf(1,'rank %d\n',j), end
    connac=zeros(m,m);
    for iloop=1:nloop;
        incr=incr+1;
        if verbose fprintf(1,' iteration %d\n',iloop), end

        % Joe: use NTU code instead (it is faster).
        % Note transposition of 'a' and 'h' since Brunet's
        % convention is that columns are samples.
        % [w,h]=nmf(a,j,verbose);
        [ h w ] = nmf_taiwan(a',rand(m,j),rand(j,n),tol,timelimit,maxiter);
        clear w; % Conserve memory.
        h = h';

        %
        % compute consensus matrix
        %
        conn=nmfconnectivity(h);
        connac=connac+conn; % accumulate connectivity matrices
    end
    consensus(j, :, :)=connac/nloop; %average
end



function conn= nmfconnectivity(h)
%
% Jean-Philippe Brunet
% Cancer Genomics 6/10/03
%
mm = size(h);
k  = mm(1);
m  = mm(2);


% compute m x m matrix which is 1 if samples are together, 0 elsewhere

% determine sample assignment by its largest metagene expresion value
[y,index] = max(h,[],1);

mat1      = repmat(index,m,1); % spread index down
mat2      = repmat(index',1,m); % spread index right

conn      = mat1==mat2; % 1 when for pair of samples with same assignement



