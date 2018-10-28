function [] = main()

% Make list of normalization methods.
%normMethods{1} = 'NoNorm';
%normMethods{2} = 'SubtractGlobalMean';
%normMethods{3} = 'ScaleColumnsThenNormalizeRows';
%normMethods{4} = 'Mean0Std1ByRows';
%normMethods{5} = 'Mean0Std1ByColumns';
%normMethods{6} = 'SubtractGlobalMean';
%normMethods{7} = 'SubtractMeanByRows';
%normMethods{8} = 'SubtractMeanByColumns';
%normMethods{9} = 'SubtractMeanByRowsAndThenByColumns';

% Make list of methods to deal with negative values.
%negMethods{1} = @SubtractAbsoluteMinimum;
%negMethods{2} = @FoldDataByRows;
%negMethods{3} = @FoldDataByColumns;
%negMethods{4} = @ExponentialScale;

% Parameters for mitEstimateK_ntuNMF_Multires.
kstart          = 2;
kend            = 20;
subsamplingRate = 4;
cophThreshold   = 0.85;
numRuns         = 50;

% Loop over normalization methods.
numNormMethods  = 9;
numNegMethods   = 4;
normEstkResults = zeros(numNormMethods,8); % Don't include Brunet, so 7 only; just 1 neg method for now.
for norm = 1:numNormMethods

    switch norm
        case 1
            normMethod = 'NoNorm';
        case 2
            normMethod = 'SubtractGlobalMean';
        case 3
            normMethod = 'ScaleColumnsThenNormalizeRows';
        case 4
            normMethod = 'Mean0Std1ByRows';
        case 5
            normMethod = 'Mean0Std1ByColumns';
        case 6
            normMethod = 'SubtractGlobalMean';
        case 7
            normMethod = 'SubtractMeanByRows';
        case 8
            normMethod = 'SubtractMeanByColumns';
        case 9
            normMethod = 'SubtractMeanByRowsAndThenByColumns';
        otherwise
            disp('Error in first switch statement')
            exit
    end

    % Just 1 neg method for now: SubtractAbsoluteMinimum
    for neg = 1:1
        % Load data.
        fileName = sprintf('/data/jmm97/NormalizationStudy/golub.%s.mat',normMethod);
        system(sprintf('touch /data/jmm97/NormalizationStudy/golub.%s.chk',normMethod));
        load(fileName);

        % Force data to be positive.
        switch neg
            case 1
                golub2 = SubtractAbsoluteMinimum(golub);
            case 2
                golub2 = FoldDataByRows(golub);
            case 3
                golub2 = FoldDataByColumns(golub);
            case 4
                golub2 = ExponentialScale(golub);
            otherwise
                disp('Error in second switch statement')
                exit
        end
        clear golub; % Conserve memory.

        % Estimate 'k' using Velicer's MAP.
        k_velicer = velicer(golub2);

        % Estimate 'k' Fogel and Young's volume-based method, with subsampling factor set to 4.
        % I.e., the interval [2,m] is divided into 4 intervals at the
        % lowest resolution level.
        kHatvec = FYandBICMultiresStrategy(golub2,4);

        k_FY    = kHatvec(1);
        k_BIC1  = kHatvec(2);
        k_BIC2  = kHatvec(3);
        k_BIC3  = kHatvec(4);
        k_rrssq = kHatvec(5);

        % Estimate 'k' using Minka's bic_pca method.
        % Commented out; this method does not appear appropriate for "wide" data.
        % [ k_minka_bic p ] = bic_pca(A);

        % Estimate 'k' using Minka's laplace_pca method.
        [ k_minka_laplace p ] = laplace_pca(golub2);

        % Try Brunet's cophenetic coefficient method.
        % Don't start at kstart=1 because the NTU NMF code doesn't like it.
        % Note transposition of data matrix 'A', because Brunet's code assumes that
        % rows are samples and columns are variables.
%        [ k_brunet cophCoefVector alreadyTried ] = mitEstimateK_ntuNMF_Multires(golub2',subsamplingRate,cophThreshold,numRuns);

        % Store results in normEstkResults
        normEstkResults(norm,1) = k_velicer;
        normEstkResults(norm,2) = k_FY;
        normEstkResults(norm,3) = k_BIC1;
        normEstkResults(norm,4) = k_BIC2;
        normEstkResults(norm,5) = k_BIC3;
        normEstkResults(norm,6) = k_rrssq;
        normEstkResults(norm,7) = k_minka_laplace;
%        normEstkResults(norm,8) = k_brunet;
    end ; % FOR neg
end ; % FOR norm

% Save results to disk.
save /data/jmm97/NormalizationStudy/normEstkResults7.mat normEstkResults

return
