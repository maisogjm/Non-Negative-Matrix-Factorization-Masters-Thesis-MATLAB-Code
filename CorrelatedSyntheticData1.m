function [ simData ] = CorrelatedSyntheticData1(m,n,k,sparsity,noiseMag,U)
% m        : Number of rows.
% n        : Number of columns.
% k        : Number of components.
% sparsity : Sparsity level (e.g., 0.4 for 40% sparsity)
% noiseMag : magnitude of the noise (e.g., 0.05 for 5% of the average magnitude of elements in simulated data
% U        : Cholesky decomposition of the Correlation Matrix C; U = chol(C);

% we randomly constructed k× matrix H with 40% sparsity. (Park and Kim)
H                  = rand(k,n); % Uniformly distributed in [0,1]
keepH              = rand(k,n); % Uniformly distributed in [0,1]
H(keepH<=sparsity) = 0.0;

% Set W to a series of orthogonal step functions.
% (Following Cichoki's "signal processing" approach.)
W = StepFunctionColumns(m,k);

% Compute the Cholesky decomposition of the correlation matrix.
% U = chol(C);

% Apply the correlations to the data to the W matrix.
W = W * U;

% Compute A (Park and Kim, and also Cichoki, construct A in this way).
A = W * H;

% and added Gaussian noise to each element where the standard
% deviation is 5% of the average magnitude of elements in A.  (Park and Kim)
avgMagA = mean(A(:));              % Compute average of elements of A.
mag     = noiseMag * avgMagA;      % 5% of the average magnitude of elements in A
A       = A + ( mag * randn(m,n)); % add Gaussian noise
simData = abs(A);                  % Make sure that A is non-negative after adding Gaussian noise.

