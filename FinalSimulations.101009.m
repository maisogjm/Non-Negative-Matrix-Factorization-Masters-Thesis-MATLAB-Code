function [] = main ()

% Obtain information from environment variables.
numSims  = str2num(getenv('NUM_SIMS'));
hostname = getenv('HOST');
run      = getenv('RUN');

% Set number of rows and columns.
m = 1000;  % Number of rows.
n = 60;    % Number of columns.

% Parameters for mitEstimateK_ntuNMF_Multires.
kstart          = 2;
kend            = 20;
subsamplingRate = 4;
cophThreshold   = 0.85;
numRuns         = 50;
%numRuns         = 1;

% Initialize output matrices.
sumsOfKhats            = zeros(1,9);
sumsOfSquaredKhats     = zeros(1,9);
sumElapsedTimes        = zeros(1,9);
sumSquaredElapsedTimes = zeros(1,9);

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

% Loop over multiple simulation runs.
for sim = 1:numSims
    % Generate simulated data using Karthik's method..
    % 1000 rows, 60 columns, 3 components
    A = SimulateData_KD();

    % Estimate 'k' using Velicer's MAP.
    tic
    k_velicer   = velicer(A);
    elapsedTime = toc;

    sumElapsedTimes(1)        = sumElapsedTimes(1)        + elapsedTime;
    sumSquaredElapsedTimes(1) = sumSquaredElapsedTimes(1) + (elapsedTime*elapsedTime);

    % Estimate 'k' Fogel and Young's volume-based method, with subsampling factor set to 4.
    % I.e., the interval [2,m] is divided into 4 intervals at the lowest resolution level.
    % Use Kullback-Liebler-based NMF computation method of Brunet et al.
    % Do numRuns re-initializations, and choose the value of k that appears most often
    % (the mode).
    rep_FY    = zeros(1,numRuns);
    rep_BIC1  = zeros(1,numRuns);
    rep_BIC2  = zeros(1,numRuns);
    rep_BIC3  = zeros(1,numRuns);
    rep_rrssq = zeros(1,numRuns);
    tic
    for rep = 1:numRuns
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.FinalSimulations.%s.%s.sim.%d.rep.%d.chk',hostname,run,sim,rep)); % check file to see progress
        kHatvec        = FYandBICMultiresStrategyKL(A,4);

        rep_FY(rep)    = kHatvec(1);
        rep_BIC3(rep)  = kHatvec(2);
        rep_BIC3(rep)  = kHatvec(3);
        rep_BIC3(rep)  = kHatvec(4);
        rep_rrssq(rep) = kHatvec(5);
    end
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.Here1.%s.%d.%d.FinalSimulations.Only1Run.chk',hostname,sim,rep)); % check file to see progress

    k_FY         = mode(rep_FY);
    k_BIC1       = mode(rep_BIC1);
    k_BIC2       = mode(rep_BIC2);
    k_BIC3       = mode(rep_BIC3);
    k_rrssq      = mode(rep_rrssq);
    elapsedTime  = toc;
    elapsedTime2 = elapsedTime*elapsedTime;
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.Here2.%s.%d.%d.FinalSimulations.Only1Run.chk',hostname,sim,rep)); % check file to see progress

    sumElapsedTimes(2)        = sumElapsedTimes(2)        + elapsedTime;
    sumElapsedTimes(3)        = sumElapsedTimes(3)        + elapsedTime;
    sumElapsedTimes(4)        = sumElapsedTimes(4)        + elapsedTime;
    sumElapsedTimes(5)        = sumElapsedTimes(5)        + elapsedTime;
    sumElapsedTimes(6)        = sumElapsedTimes(6)        + elapsedTime;
    sumSquaredElapsedTimes(2) = sumSquaredElapsedTimes(2) + elapsedTime2;
    sumSquaredElapsedTimes(3) = sumSquaredElapsedTimes(3) + elapsedTime2;
    sumSquaredElapsedTimes(4) = sumSquaredElapsedTimes(4) + elapsedTime2;
    sumSquaredElapsedTimes(5) = sumSquaredElapsedTimes(5) + elapsedTime2;
    sumSquaredElapsedTimes(6) = sumSquaredElapsedTimes(6) + elapsedTime2;

    % Estimate 'k' using Minka's bic_pca method.
    % Commented out; this method does not appear appropriate for "wide" data.
    % [ k_minka_bic p ] = bic_pca(A);

%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.Here3.%s.%d.%d.FinalSimulations.Only1Run.chk',hostname,sim,rep)); % check file to see progress
    % Estimate 'k' using Minka's laplace_pca method.
    tic
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.laplace_pca.%s.FinalSimulations.Only1Run.chk',hostname)); % check file to see progress
    [ k_minka_laplace p ] = laplace_pca(A);
    elapsedTime           = toc;

%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.Here3.%s.%d.%d.FinalSimulations.Only1Run.chk',hostname,sim,rep)); % check file to see progress
    sumElapsedTimes(7)        = sumElapsedTimes(7)        + elapsedTime;
    sumSquaredElapsedTimes(7) = sumSquaredElapsedTimes(7) + (elapsedTime*elapsedTime);

    % Try Brunet's cophenetic coefficient method.
    % Don't start at kstart=1 because the NTU NMF code doesn't like it.
    % Note transposition of data matrix 'A', because Brunet's code assumes that
    % rows are samples and columns are variables.
    % Use Kullback-Liebler-based NMF computation method of Brunet et al.
    tic
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.CCC.%s.FinalSimulations.Only1Run.chk',hostname)); % check file to see progress
    [ k_brunet cophCoefVector alreadyTried ] = mitEstimateK_mitNMF_Multires2(A',subsamplingRate,cophThreshold,numRuns);
    elapsedTime                              = toc;

    sumElapsedTimes(8)        = sumElapsedTimes(8)        + elapsedTime;
    sumSquaredElapsedTimes(8) = sumSquaredElapsedTimes(8) + (elapsedTime*elapsedTime);

    % Try Owen and Perry's BCV method.
    tic
    kList = [ 2:10 ];
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.BCV1.%s.FinalSimulations.Only1Run.chk',hostname)); % check file to see progress
    [ rowgps, colgps ] = genbicvgps(m, n, 2, 5);
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.BCV2.%s.FinalSimulations.Only1Run.chk',hostname)); % check file to see progress
    [ gI ]             = gengabcvsvd( A,  [ rowgps, colgps ], kList, 0);
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.BCV3.%s.FinalSimulations.Only1Run.chk',hostname)); % check file to see progress
    elapsedTime        = toc;
    [ val valIndex ]   = max(gI);
    k_BCV              = kList(valIndex);

    sumElapsedTimes(9)        = sumElapsedTimes(9)        + elapsedTime;
    sumSquaredElapsedTimes(9) = sumSquaredElapsedTimes(9) + (elapsedTime*elapsedTime);

    % Update sums and sums of squares.
    kHatVector                   = [ k_velicer k_FY k_BIC1 k_BIC1 k_BIC1 k_rrssq k_minka_laplace k_brunet k_BCV ];
    sumsOfKhats        = sumsOfKhats + kHatVector
    sumsOfSquaredKhats = sumsOfKhats + (kHatVector.*kHatVector);
end ; % FOR sim
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.FinishedLoop.%s.FinalSimulations.Only1Run.chk',hostname)); % check file to see progress

% Save results to disk.
% Include the random number seed that was used to generate the simulation data.
% Include the total time for the MATLAB session.
totalTime = cputime;
outfile   = getenv('OUTFILE');
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.SaveTo.%s.%s.FinalSimulations.Only1Run.chk',outfile,hostname)); % check file to see progress
save(outfile,'seed','numSims','totalTime','sumsOfKhats','sumsOfSquaredKhats','sumElapsedTimes','sumSquaredElapsedTimes')

% Save the current state of the random number seed to disk.
% This will overwrite the previous file on disk.
seed = rand('state',seed); 
save(seedFile,'seed')

% Exit MATLAB.
exit

return
