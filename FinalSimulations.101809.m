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
kList           = [ kstart:kend ];
subsamplingRate = 4;
cophThreshold   = 0.85;
numRuns         = 50;
%numRuns         = 1;

% Initialize output matrices.
k_velicer       = zeros(numSims,1); % Non-iterative
k_minka_laplace = zeros(numSims,1); % Non-iterative
k_minka_bic     = zeros(numSims,1); % Non-iterative
k_FY            = zeros(numSims,numRuns);
k_BIC1          = zeros(numSims,numRuns);
k_BIC2          = zeros(numSims,numRuns);
k_BIC3          = zeros(numSims,numRuns);
k_rrssq         = zeros(numSims,numRuns);
k_CCC           = zeros(numSims,kend-kstart+1);
k_BCV           = zeros(numSims,kend-kstart+1);

t_velicer       = zeros(numSims,1); % Non-iterative
t_minka_laplace = zeros(numSims,1); % Non-iterative
t_minka_bic     = zeros(numSims,1); % Non-iterative
t_FY            = zeros(numSims,numRuns);
t_BIC1          = zeros(numSims,numRuns);
t_BIC2          = zeros(numSims,numRuns);
t_BIC3          = zeros(numSims,numRuns);
t_rrssq         = zeros(numSims,numRuns);
t_CCC           = zeros(numSims,kend-kstart+1);
t_BCV           = zeros(numSims,kend-kstart+1);

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
    A = SimulateData_KD3();

    % Estimate 'k' using Velicer's MAP.
    tic
    k_velicer(sim) = velicer(A);
    t_velicer(sim) = toc;

    % Estimate 'k' Fogel and Young's volume-based method, with subsampling factor set to 4.
    % I.e., the interval [2,m] is divided into 4 intervals at the lowest resolution level.
    % Use Kullback-Liebler-based NMF computation method of Brunet et al.
    % Do numRuns re-initializations, and choose the value of k that appears most often
    % (the mode).
    for rep = 1:numRuns
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.FinalSimulations.%s.%s.sim.%d.rep.%d.chk',hostname,run,sim,rep)); % check file to see progress
        tic
        kHatvec          = FYandBICMultiresStrategyKL(A,4);
        elapsed_time     = toc

        k_FY(sim,rep)    = kHatvec(1);
        k_BIC3(sim,rep)  = kHatvec(2);
        k_BIC3(sim,rep)  = kHatvec(3);
        k_BIC3(sim,rep)  = kHatvec(4);
        k_rrssq(sim,rep) = kHatvec(5);

        t_FY(sim,rep)    = elapsed_time;
        t_BIC3(sim,rep)  = elapsed_time;
        t_BIC3(sim,rep)  = elapsed_time;
        t_BIC3(sim,rep)  = elapsed_time;
        t_rrssq(sim,rep) = elapsed_time;
    end
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.Here1.%s.%d.%d.FinalSimulations.Only1Run.chk',hostname,sim,rep)); % check file to see progress

    % Estimate 'k' using Minka's laplace_pca method.
    tic
    [ k_minka_laplace(sim) p ] = laplace_pca(A');
    t_minka_laplace(sim)       = toc;

    % Estimate 'k' using Minka's bic_pca method.
    tic
    [ k_minka_bic(sim) p ] = bic_pca(A);
    t_minka_bic(sim)       = toc;

    % Try Brunet's cophenetic coefficient method.
    % Don't start at kstart=1 because the NTU NMF code doesn't like it.
    % No longer using NTU NMF code, but don't start at kstart=1 anyway just in case.
    % Use Kullback-Liebler-based NMF computation method of Brunet et al.
    % Note transposition of data matrix 'A', because Brunet's code assumes that
    % rows are samples and columns are variables.
    % Actually, this transposition may be unnecessary.
    tic
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.CCC.%s.FinalSimulations.Only1Run.chk',hostname)); % check file to see progress
%    [ k_brunet cophCoefVector alreadyTried ] = mitEstimateK_mitNMF_Multires2(A',subsamplingRate,cophThreshold,numRuns);
    consensus                            = nmfconsensus(A',kstart,kend,numRuns,0);
    [ ordcons, clustid, ordindex, coph ] = nmforderconsensus(consensus,kstart,kend);
    elapsedTime                          = toc;
    k_CCC(sim,:) = coph(:);

    % Try Owen and Perry's BCV method.
    tic
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.BCV1.%s.FinalSimulations.Only1Run.chk',hostname)); % check file to see progress
    [ rowgps, colgps ] = genbicvgps(m, n, 2, 5);
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.BCV2.%s.FinalSimulations.Only1Run.chk',hostname)); % check file to see progress
    [ gI ]             = gengabcvsvd( A,  [ rowgps, colgps ], kList, 0);
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.BCV3.%s.FinalSimulations.Only1Run.chk',hostname)); % check file to see progress
    elapsedTime        = toc;
    k_BCV(sim,:)       = gI(:);

end ; % FOR sim

%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.FinishedLoop.%s.FinalSimulations.Only1Run.chk',hostname)); % check file to see progress

% Save results to disk.
% Include the random number seed that was used to generate the simulation data.
% Include the total time for the MATLAB session.
totalTime = cputime;
outfile   = getenv('OUTFILE');
%system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/Sim.SaveTo.%s.%s.FinalSimulations.Only1Run.chk',outfile,hostname)); % check file to see progress
save(outfile,'seed','numSims','kList','totalTime', ...
    'k_velicer','k_minka_laplace','k_minka_bic','k_FY','k_BIC1','k_BIC2','k_BIC3','k_rrssq','k_CCC','k_BCV', ...
    't_velicer','t_minka_laplace','t_minka_bic','t_FY','t_BIC1','t_BIC2','t_BIC3','t_rrssq','t_CCC','t_BCV')

% Save the current state of the random number seed to disk.
% This will overwrite the previous file on disk.
seed = rand('state',seed); 
save(seedFile,'seed')

% Exit MATLAB.
exit

return
