function [GoodClustersData, Mean_GB,STD_GB,nb_GB] = Clusters_Data(Goodbetas,idxKmeans_DF,DF)

for i=1:length(Goodbetas)    
    GoodClustersData(i).DF=DF(idxKmeans_DF==Goodbetas(i),:)*100;
    GoodClustersData(i).Mean=mean(GoodClustersData(i).DF,1);
    GoodClustersData(i).STD=std(GoodClustersData(i).DF,1,1);    
end

Mean_GB=[GoodClustersData.Mean];Mean_GB=reshape(Mean_GB,length(GoodClustersData(1).Mean),length(GoodClustersData));
STD_GB=[GoodClustersData.STD];STD_GB=reshape(STD_GB,length(GoodClustersData(1).STD),length(GoodClustersData));
nb_GB=ones(length(GoodClustersData),length(GoodClustersData(1).Mean));
for i=1:length(Goodbetas)
    nb_GB(i,:)=size(GoodClustersData(i).DF,1)*nb_GB(i,:);
end
nb_GB=nb_GB';