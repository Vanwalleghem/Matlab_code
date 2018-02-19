DF=DeltaF2(double(SelectCorr),21,11);
max_DF=max(DF,[],2);
min_DF=min(DF,[],2);
idx_Threshold_5to200=find(max_DF>0.05 & max_DF<2 & min_DF>-0.1);
idx_bool_Thr=(max_DF>0.05 & max_DF<2 & min_DF>-0.1);
Select_DF_Thr5=DF(idx_Threshold_5to200,:);

ZS=zscore(Calcium,1,2);
random_5k=datasample(ZS,5000);
eva = evalclusters(random_5k,'kmeans','silhouette','Distance','correlation','KList',[1:100]);
options = statset('UseParallel',1); [idxKmeans Cmap]=kmeans(ZS,eva.OptimalK,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');

