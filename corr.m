function [ C ] = corr(I)

%[ N M ] = size(I);

%globalMean = mean(I,1);
%for i=1:M
%    I(:,i) = I(:,i)-globalMean(i);
%end
%standardDev = std(I,1);
%for i=1:M
%    I(:,i) = I(:,i)/standardDev(i);
%end

%C = I'*I;
C = corrcoef(I);
