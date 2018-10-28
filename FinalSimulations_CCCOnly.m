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

% Try Brunet's cophenetic coefficient method.
% Use Kullback-Liebler-based NMF computation method of Brunet et al.
% Note transposition of data matrix 'A', because Brunet's code assumes that
% rows are samples and columns are variables.
% This transposition may be unnecessary.
system(sprintf('touch /data/jmm97/PrototypeBigSimulation3/progress/FinalSimulations/Sim.%s.%s.%s.%d.CCC.chk',vers,hostname,run,trueK)); % check file to see progress
tic
consensus                            = nmfconsensus(A',kstart,kend,numRuns,0);
[ ordcons, clustid, ordindex, coph ] = nmforderconsensus(consensus,kstart,kend);
elapsedTime                          = toc;
k_CCC = coph(:);
t_CCC = elapsedTime;

% Save results to disk.
% Include the random number seed that was used to generate the simulation data.
% Include the total time for the MATLAB session.
totalTime = cputime;
save(outfile,'seed','numRuns','kList','trueK','totalTime', 'k_CCC','t_CCC')

% Save the current state of the random number seed to disk.
% This will overwrite the previous file on disk.
rand('state',seed); 
save(seedFile,'seed')

% Exit MATLAB.
exit

return
