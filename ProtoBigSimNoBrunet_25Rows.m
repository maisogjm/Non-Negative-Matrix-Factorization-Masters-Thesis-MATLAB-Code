function [] = main ()

% Obtain information from environment variables.
numSims  = str2num(getenv('NUM_SIMS'));
hostname = getenv('HOST');
runLabel = getenv('RUN');

% We used m = 25 and n = 1000.
m = 25;   % Number of rows.
n = 1000; % Number of columns.

% Parameters for mitEstimateK_ntuNMF_Multires.
kstart          = 2;
kend            = 20;
subsamplingRate = 4;
cophThreshold   = 0.85;
numRuns         = 50;

% Define list of experimental k's.
%kList = [ 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 49, 50, 60, 70, 80, 88, 89, 90, 100 ];
kList = 2:20;
numK  = length(kList);

% Initialize output matrices.
sumsOfKhats            = zeros(numK,7);
sumsOfSquaredKhats     = zeros(numK,7);
sumElapsedTimes        = zeros(numK,7);
sumSquaredElapsedTimes = zeros(numK,7);

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
%disp(sprintf('k = %d\n',k))
[s,w]  = system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.%s.%d.%s.chk',hostname,kIndex,runLabel)); % check file to see progress
    for sim = 1:numSims
        % Generate simulated data.
        % m rows, n columns, k components
        % 40% sparsity in the H matrix, 5% Gaussian noise
        A = ParkKimCichokiHybridSimData(m,n,kIndex,0.4,0.05);

        % Estimate 'k' using Velicer's MAP.
        tic
        k_velicer   = velicer(A);
        elapsedTime = toc;

        sumElapsedTimes(kIndex,1)        = sumElapsedTimes(kIndex,1)        + elapsedTime;
        sumSquaredElapsedTimes(kIndex,1) = sumSquaredElapsedTimes(kIndex,1) + (elapsedTime*elapsedTime);

        % Estimate 'k' Fogel and Young's volume-based method, with subsampling factor set to 4.
        % I.e., the interval [2,m] is divided into 4 intervals at the lowest resolution level.
%[s,w]  = system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.%s.%d.%d.%s.FY_BIC.chk',hostname,kIndex,sim,runLabel)); % check file to see progress
        tic
        kHatvec      = FYandBICMultiresStrategy(A,4);
        elapsedTime  = toc;
        elapsedTime2 = elapsedTime*elapsedTime;

        k_FY         = kHatvec(1);
        k_BIC1       = kHatvec(2);
        k_BIC2       = kHatvec(3);
        k_BIC3       = kHatvec(4);
        k_rrssq      = kHatvec(5);

        sumElapsedTimes(kIndex,2)        = sumElapsedTimes(kIndex,2)        + elapsedTime;
        sumElapsedTimes(kIndex,3)        = sumElapsedTimes(kIndex,3)        + elapsedTime;
        sumElapsedTimes(kIndex,4)        = sumElapsedTimes(kIndex,4)        + elapsedTime;
        sumElapsedTimes(kIndex,5)        = sumElapsedTimes(kIndex,5)        + elapsedTime;
        sumElapsedTimes(kIndex,6)        = sumElapsedTimes(kIndex,6)        + elapsedTime;
        sumSquaredElapsedTimes(kIndex,2) = sumSquaredElapsedTimes(kIndex,2) + elapsedTime2;
        sumSquaredElapsedTimes(kIndex,3) = sumSquaredElapsedTimes(kIndex,3) + elapsedTime2;
        sumSquaredElapsedTimes(kIndex,4) = sumSquaredElapsedTimes(kIndex,4) + elapsedTime2;
        sumSquaredElapsedTimes(kIndex,5) = sumSquaredElapsedTimes(kIndex,5) + elapsedTime2;
        sumSquaredElapsedTimes(kIndex,6) = sumSquaredElapsedTimes(kIndex,6) + elapsedTime2;

        % Estimate 'k' using Minka's bic_pca method.
        % Commented out; this method does not appear appropriate for "wide" data.
        % [ k_minka_bic p ] = bic_pca(A);

        % Estimate 'k' using Minka's laplace_pca method.
        tic
        [ k_minka_laplace p ] = laplace_pca(A);
        elapsedTime           = toc;

        sumElapsedTimes(kIndex,7)        = sumElapsedTimes(kIndex,7)        + elapsedTime;
        sumSquaredElapsedTimes(kIndex,7) = sumSquaredElapsedTimes(kIndex,7) + (elapsedTime*elapsedTime);

        % Try Brunet's cophenetic coefficient method.
        % Don't start at kstart=1 because the NTU NMF code doesn't like it.
        % Note transposition of data matrix 'A', because Brunet's code assumes that
        % rows are samples and columns are variables.
%[s,w]  = system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.%d.%d.Brunet.chk',kIndex,sim)); % check file to see progress
%        tic
%        [ k_brunet cophCoefVector alreadyTried ] = mitEstimateK_ntuNMF_Multires(A',subsamplingRate,cophThreshold,numRuns);
%        elapsedTime                              = toc;

%        sumElapsedTimes(kIndex,8)        = sumElapsedTimes(kIndex,8)        + elapsedTime;
%        sumSquaredElapsedTimes(kIndex,8) = sumSquaredElapsedTimes(kIndex,8) + (elapsedTime*elapsedTime);

        % Update sums and sums of squares.
        kHatVector                   = [ k_velicer k_FY k_BIC1 k_BIC1 k_BIC1 k_rrssq k_minka_laplace ];
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
save(outfile,'seed','numSims','kList','totalTime','sumsOfKhats','sumsOfSquaredKhats','sumElapsedTimes','sumSquaredElapsedTimes')

% Save the current state of the random number seed to disk.
% This will overwrite the previous file on disk.
rand('state',seed); 
save(seedFile,'seed')

% Exit MATLAB.
exit

return
