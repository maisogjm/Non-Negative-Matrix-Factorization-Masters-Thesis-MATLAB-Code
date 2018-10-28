function [] = main ()

%[s,w]  = system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/EnteredFunction.chk')); % check file to see progress
% Obtain information from environment variables.
numSims = str2num(getenv('NUM_SIMS'));
hostname = getenv('HOST');
runLabel = getenv('RUN');

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
% kList = [ 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 49, 50, 60, 70, 80, 88, 89, 90, 100 ];
kList = 2:20;
numK  = length(kList);

%[s,w]  = system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/GotThisFar1.chk')); % check file to see progress
% Initialize output matrices.
sumsOfKhats            = zeros(numK,8);
sumsOfSquaredKhats     = zeros(numK,8);
sumElapsedTimes        = zeros(numK,8);
sumSquaredElapsedTimes = zeros(numK,8);
cophCoefVectorMAT      = zeros(numK,numSims);

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
kIndex = 1;
for k = kList
[s,w]  = system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/NewBrunetOnly.%s.%d.%s.chk',hostname,k,runLabel)); % check file to see progress
    % Set up correlation matrix, compute its Cholesky decomposition U.
    offDiagVals = 0.5;
    C           = GenerateCorrelationMatrix(k,offDiagVals);
    U           = chol(C);

    for sim = 1:numSims
        % Generate simulated data.
        % m rows, n columns, k components
        % 40% sparsity in the H matrix, 5% Gaussian noise
        % A = ParkKimCichokiHybridSimData(m,n,k,0.4,0.05);
        A = CorrelatedSyntheticData1(m,n,k,0.4,0.05,U);

        % Try Brunet's cophenetic coefficient method.
        % Don't start at kstart=1 because the NTU NMF code doesn't like it.
        tic
        [ k_brunet cophCoefVector alreadyTried ] = mitEstimateK_ntuNMF_Multires2(A,subsamplingRate,cophThreshold,numRuns);
        elapsedTime                              = toc;

        sumElapsedTimes(kIndex,8)        = sumElapsedTimes(kIndex,8)        + elapsedTime;
        sumSquaredElapsedTimes(kIndex,8) = sumSquaredElapsedTimes(kIndex,8) + (elapsedTime*elapsedTime);
        alreadyTriedTMP                  = find(alreadyTried~=0);
        cophCoefVectorMAT(:,sim)         = cophCoefVector(alreadyTriedTMP);

        % Update sums and sums of squares.
        kHatVector                   = [ 0, 0, 0, 0, 0, 0, 0, k_brunet ];
save /data/jmm97/PrototypeBigSimulation3/progress/temp.mat
        sumsOfKhats(kIndex,:)        = sumsOfKhats(kIndex,:) + kHatVector;
        sumsOfSquaredKhats(kIndex,:) = sumsOfSquaredKhats(kIndex,:) + (kHatVector.*kHatVector);
    end ; % FOR sim

    kIndex = kIndex + 1;
end ; % FOR k

% Save results to disk.
% Include the random number seed that was used to generate the simulation data.
% Include the total time for the MATLAB session.
totalTime = cputime;
outfile   = getenv('OUTFILE');
save(outfile,'seed','numSims','kList','totalTime','sumsOfKhats','sumsOfSquaredKhats','sumElapsedTimes','sumSquaredElapsedTimes','cophCoefVectorMAT')

% Save the current state of the random number seed to disk.
% This will overwrite the previous file on disk.
rand('state',seed); 
save(seedFile,'seed')

% Exit MATLAB.
exit

return
