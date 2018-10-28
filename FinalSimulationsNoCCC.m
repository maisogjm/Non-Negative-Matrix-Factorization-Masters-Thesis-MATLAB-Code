function [] = main ()

% Obtain information from environment variables.
hostname = getenv('HOST');
run      = getenv('RUN');
outfile  = getenv('OUTFILE');
seedFile = getenv('SEEDFILE');
vers     = getenv('VERSION'); % 'version' is a built-in MATLAB command
trueK    = str2num(getenv('NUM_K'));

% Set number of rows and columns.
m = 1000;  % Number of rows.
n = 60;    % Number of columns.

% Parameters for mitEstimateK_ntuNMF_Multires.
kstart          = 2;
kend            = 20;
kList           = [ kstart:kend ];
cophThreshold   = 0.85;
numRuns         = 50;
%numRuns         = 1;

% Initialize output matrices.
k_velicer       = zeros(numRuns,1); % Non-iterative
k_minka_laplace = zeros(numRuns,1); % Non-iterative
k_minka_bic     = zeros(numRuns,1); % Non-iterative
k_FY            = zeros(numRuns,kend-kstart+1);
k_BIC1          = zeros(numRuns,kend-kstart+1);
k_BIC2          = zeros(numRuns,kend-kstart+1);
k_BIC3          = zeros(numRuns,kend-kstart+1);
k_rrssq         = zeros(numRuns,kend-kstart+1);
k_BCV           = zeros(numRuns,kend-kstart+1);

t_velicer       = zeros(numRuns,1); % Non-iterative
t_minka_laplace = zeros(numRuns,1); % Non-iterative
t_minka_bic     = zeros(numRuns,1); % Non-iterative
t_FY            = zeros(numRuns,1);
t_BIC1          = zeros(numRuns,1);
t_BIC2          = zeros(numRuns,1);
t_BIC3          = zeros(numRuns,1);
t_rrssq         = zeros(numRuns,1);
t_BCV           = zeros(numRuns,1);

% Initialize random number seed.
% Attempt to load the random number seed from a seed file on disk.
try
    load(seedFile)

% If unsuccessful, generate the seed from the clock time.
% The idea of summing 100 times the clock is from the MATLAB documentation for RAND.
catch
    seed = sum(100*clock);
end
rand('state',seed); % Set the random number seed.

% Generate simulated data using Karthik's method.
% 1000 rows, 60 columns.  # of components = trueK.
A = SimulateData_KD_N(trueK);

% Loop over multiple simulation runs.
for rep = 1:numRuns

system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/Sim.%s.%s.%d.%d.chk',vers,hostname,run,trueK,rep)); % check file to see progress

    % Estimate 'k' using Velicer's MAP.
    tic
    k_velicer(rep) = velicer(A);
    t_velicer(rep) = toc;

    % Estimate 'k' using Minka's laplace_pca method.
    % Seems to work better with the matrix transposed.
    tic
    [ k_minka_laplace(rep) p ] = laplace_pca(A');
    t_minka_laplace(rep)       = toc;

    % Estimate 'k' using Minka's bic_pca method.
    % Doesn't seem to work so well, but let's compute it anyway.
    tic
    [ k_minka_bic(rep) p ] = bic_pca(A);
    t_minka_bic(rep)       = toc;

    % Estimate 'k' Fogel and Young's volume-based method
    % and the three BIC methods.  No subsampling.
    % Use Kullback-Liebler-based NMF computation method of Brunet et al.
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/InvokingFY_BIC_RRSSQ_KL.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));
    tic
    [ kVolume kBIC1 kBIC2 kBIC3 kRRSSQ ] = FY_BIC_RRSSQ_KL(A,kList);
    elapsed_time   = toc;

system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/Volume.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));
    k_FY(rep,:)    = kVolume(:)';
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/BIC1.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));
    k_BIC1(rep,:)  = kBIC1(:)';
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/BIC2.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));
    k_BIC2(rep,:)  = kBIC2(:)';
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/BIC3.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));
    k_BIC3(rep,:)  = kBIC3(:)';
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/RRSSQ.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));
    k_rrssq(rep,:) = kRRSSQ(:)';
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/ALLDONE.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));

    t_FY(rep)    = elapsed_time;
    t_BIC1(rep)  = elapsed_time;
    t_BIC2(rep)  = elapsed_time;
    t_BIC3(rep)  = elapsed_time;
    t_rrssq(rep) = elapsed_time;

    % Try Owen and Perry's BCV method.
    % 2-fold cross-validation over the 1000 rows.
    % 5-fold cross-validation over the 60 columns.
    tic
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/genbicvgps.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));
    [ rowgps, colgps ] = genbicvgps(m, n, 2, 5);
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/gengabcvsvd.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));
    [ gI ]             = gengabcvsvd( A,  rowgps, colgps, kList, 0);
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/toc.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));
    elapsedTime        = toc;
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/gI.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));
    k_BCV(rep,:)       = gI(:);
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/elapsedTime.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));
    t_BCV(rep)         = elapsedTime;
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/ALLDONE2.%s.%s.%d.%d.%dchk',vers,hostname,run,trueK,rep));

end ; % FOR rep

% Save results to disk.
% Include the random number seed that was used to generate the simulation data.
% Include the total time for the MATLAB session.
totalTime = cputime;
save(outfile,'seed','numRuns','kList','trueK','totalTime', ...
    'k_velicer','k_minka_laplace','k_minka_bic','k_FY','k_BIC1','k_BIC2','k_BIC3','k_rrssq','k_BCV', ...
    't_velicer','t_minka_laplace','t_minka_bic','t_FY','t_BIC1','t_BIC2','t_BIC3','t_rrssq','t_BCV')

% Save the current state of the random number seed to disk.
% This will overwrite the previous file on disk.
rand('state',seed); 
save(seedFile,'seed')

% Exit MATLAB.
exit

return
