function [] = main()

hostname = getenv('HOST');

% We used m = 100 and n = 1000.
m = 100;  % Number of rows.
n = 1000; % Number of columns.

% Parameters for mitEstimateK_ntuNMF_Multires.
kstart          = 2;
kend            = 20;
subsamplingRate = 4;
cophThreshold   = 0.99;
numRuns         = 20;

% Define list of experimental k's.
kList = 2:20;
numK  = length(kList);

% Initialize output matrices.
cophCoefVectorMAT = zeros(numK,100);
alreadyTriedMAT   = zeros(numK,100);
kHatMAT           = zeros(numK,1);

% Initialize random number seed.
% Attempt to load the random number seed from a seed file on disk.
seedFile = getenv('SEEDFILE');
try
    load(seedFile)

% If unsuccessful, generate the seed from the clock time.
% The idea of summing 100 times the clock is from the MATLAB documentation for RAND.
catch
    seed = sum(100*clock);
end
rand('state',seed); % Set the random number seed.

% Simulate different numbers of underlying components, 'k'.
%system('rm /tmp/*.chk')
for k = kList
%    disp(sprintf('k = %d',k));
%    system(sprintf('touch /tmp/BrunetOnly.%s.k.%d.CccThr.%g.chk\n',hostname,k,cophThreshold));
    % Generate simulated data.
    % m rows, n columns, k components
    % 40% sparsity in the H matrix, 5% Gaussian noise
    A = ParkKimCichokiHybridSimData(m,n,k,0.4,0.05);

    % Try Brunet's cophenetic coefficient method.
    % Don't start at kstart=1 because the NTU NMF code doesn't like it.
    [ cophCoefVector alreadyTried kHat ]  ...
   = mitEstimateK_ntuNMF_Multires2(A,subsamplingRate,cophThreshold,numRuns,k);

    cophCoefVectorMAT(k,:) = cophCoefVector(:);
    alreadyTriedMAT(k,:)   = alreadyTried(:);
    khatMAT(k)             = kHat;
end ; % FOR k

% Save results to disk.
% Include the random number seed that was used to generate the simulation data.
% Include the total time for the MATLAB session.
outFile = getenv('OUTFILE');
save(outFile,'cophCoefVectorMAT','alreadyTriedMAT','khatMAT');

% Save the current state of the random number seed to disk.
% This will overwrite the previous file on disk.
rand('state',seed); 
save(seedFile,'seed')

% Exit MATLAB.
exit

return

