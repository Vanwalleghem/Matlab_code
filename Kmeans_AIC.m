function [AIC,BIC] = Kmeans_AIC(ClustCenters,idxKmeans,sumD)
%Computes AIC and BIC of Kmeans result, requires the
%[idxKmeans,ClustCenters,sumD]=kmeans output
[k,m]=size(ClustCenters);
n=length(idxKmeans);
D=sum(sumD);
AIC=D+2*m*k;
BIC=D+log(n)*m*k;

