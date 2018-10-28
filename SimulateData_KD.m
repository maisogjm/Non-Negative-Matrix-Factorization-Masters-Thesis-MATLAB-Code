function [ simdat ] = SimulateData_KD()
%
% FORMAT [ simdat ] = SimulateData_KD()
%
% Joe Maisog's MATLAB translation of Karthik Devarajan's
% R code for creating a 1000X60 matrix of randomly generated gene expression profiles.
%
% Data consists of k=3 clusters: 
% Cluster 1 is based on the expression profiles of the first 50 genes (rows 1-50)
% across the first 20 samples (columns 1-20) generated from an exponential distribution
% with mean 40 (simdat2 below).
%
% Cluster 2 is based on the expression profiles of the first 100 genes (rows 1-100)
% across the next 20 samples (columns 21-40) generated from an exponential distribution
% with mean 80 (simdat31, simdat3).
% 
% Cluster 3 is based on the expression profiles of the second set of 50 genes (rows 51-100)
% across the next 20 samples (columns 41-60) generated from an exponential distribution
% with mean 120 (simdat4). 
%
% The remaining 900 rows represent the expression profiles of "noisy" genes generated from
% the unit exponential distribution (simdat1)

% Generate the 900 rows of "noisy" genes.       % Below are line numbers and statements in the original R code.
simdat1 = exprnd(900,60,1);                     % 17: simdat1[i,] <- rexp(60,1)

% Fill cluster submatrices with data.
simdat2   = exprnd(50,20, 1/40);                % 24: simdat2[i,] <- rexp(20,1/40)
simdat4   = exprnd(50,20, 1/120);               % 25: simdat4[i,] <- rexp(20,1/120)
simdat3   = exprnd(50,20, 1/80);                % 26: simdat3[i,] <- rexp(20,1/80)
simdat31  = exprnd(50,20, 1/80);                % 27: simdat31[i,] <- rexp(20,1/80)
simdat11  = exprnd(50,20, 1);                   % 28: simdat11[i,] <- rexp(20,1)
simdat111 = exprnd(50,20, 1);                   % 29: simdat111[i,] <- rexp(20,1)

% Construct output matrix from submatrices.
% 31: simdat <- rbind(cbind(simdat2,simdat31,simdat111),cbind(simdat11,simdat3,simdat4),simdat1)
simdat = [ [simdat2 , simdat31 , simdat111] ; [ simdat11 , simdat3 , simdat4 ] ; simdat1 ];

% writing data onto a file
% write.table(simdat,file="simdat.txt",quote=FALSE,sep="\t",col.names=FALSE,row.names=FALSE)
