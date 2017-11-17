MatFiles=dir('G*.mat');
name=strcat(MatFiles(1).name);
% Baseline=load(name, 'Baseline');
% Baseline=Baseline.Baseline;
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
Spikes=load(name, 'Spikes');
Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
Fitness=load(name, 'fitness');
Fitness=Fitness.fitness;
for i = 2:length(MatFiles)
name=strcat(MatFiles(i).name);
% B=load(name, 'Baseline');
% B=B.Baseline;
C=load(name, 'DenoisedTraces');
C=C.DenoisedTraces;
S=load(name, 'Spikes');
S=S.Spikes;
N=load(name, 'Noise');
N=N.Noise;
F=load(name, 'fitness');
F=F.fitness;
Noise=vertcat(Noise,N);
Calcium=vertcat(Calcium,C);
Spikes=vertcat(Spikes,S);
Fitness=horzcat(Fitness,F);
% Baseline=horzcat(Baseline,B);
MatFiles(i).number=size(Calcium,1);
end
clearvars C S F N name i B;

Fit_Calcium=Calcium(Fitness<-40,:);

parfor i=1:size(Select_DF_Thr5,1)
    mdl=stepwiselm(aud8freq',Select_DF_Thr5(i,:),'linear','Criterion','adjrsquared','Upper','interactions','verbose',0);
    model_DF_Thr5(i).coef=mdl.Coefficients;
    model_DF_Thr5(i).MSE=mdl.MSE;
    model_DF_Thr5(i).Fitted=mdl.Fitted;
    model_DF_Thr5(i).rsquared=mdl.Rsquared.Adjusted;
end

Coeff=zeros(length(model),8);
for i=1:length(A)
    B=A{i};B=double(B);B(B==0)=NaN;
    image=zeros(size(B,1),size(B,2),3);        
    for j=min(min(B)):max(max(B))
        temp=zeros(size(B));      
        for k=7:8
            if Model_all(j).coef.pValue(k) < 0.01              
                temp(B==j)=Model_all(j).coef.Estimate(k)*Model_all(j).rsquared;
                image(:,:,k-6)=image(:,:,k-6)+temp;
            end
        end
    end
    image=image/max(max(max(image)));image=image*256;
    image=uint8(image); 
    name=strcat('rsqColor3_',sorted_files{i});
    imwrite(image,name,'tif');
end
clearvars image i j name B;

parfor i=1:size(Cmap_DF_20,1)
    mdl=stepwiselm(aud5vol',Cmap_DF_20(i,:),'linear','Criterion','adjrsquared','Upper','interactions','verbose',0);
    model_sp(i).coef=mdl.Coefficients;
    model_sp(i).MSE=mdl.MSE;
    model_sp(i).Fitted=mdl.Fitted;
    model_sp(i).rsquared=mdl.Rsquared.Adjusted;
end
clearvars mdl i;

Numbers=[0 [MatFiles.number]];temp=find(idxKmeans_ZS==SelectSClusters(1));
AVG_files=dir('D:\Pictures\processed\Tonotropy\Tectum\AVG*.tif');
fileNb=find(Numbers<temp(1),1,'last');filename=MatFiles(fileNb).name;imagename=regexp(filename,'_output_analysis_matlab.mat','split');imagename=strcat('AVG_',imagename{1},'.tif');
image=imread(strcat('D:\Pictures\processed\Tonotropy\Tectum\',imagename));image=repmat(double(image),1,1,3);
for idx=SelectSClusters
    temp=find(idxKmeans_ZS==idx);
    for i=temp'
        tempFileNb=find(Numbers<i,1,'last');
        if tempFileNb ~= fileNb
            
        end    
    end
end
clearvars idx i temp tempFileNb fileNb AVG_files filename image

idx_final_modelKmeans=zeros(length(ZS),1);
counter=1;
for i=idxKmeans_selectZS'
    idx_final_modelKmeans(Select_model_idx(counter))=i;
    counter=counter+1;
end
clearvars i counter

Select_spikes=Spikes(~low_Std,:);
idx_final=zeros(length(Spikes),1);
counter=1;
for i=1:length(low_Std)
    if ~low_Std(i)
        idx_final(i)=idxKmeans_Spikes_20(counter);
        counter=counter+1;
    end
end
clearvars counter;

%MatFiles=dir('G*.mat');
Numbers=[0 [MatFiles.number]];
temp=[];
counter=1;
for idx=Select_model_Clusters
    temp{counter}=find(idxKmeans_Spikes_20==idx);
    %tempidx=find(idxKmeans==idx);
    %temp{counter}=Select(tempidx)';
    counter=counter+1;    
end

Start=min(cellfun(@min, temp));Start=find(Numbers<Start,1,'last');
filename=MatFiles(Start).name;

for idx=Start:length(MatFiles)
    filename=MatFiles(idx).name;
    ROIsNb=[];ClusterNb=[];
    %for k = 1 : length(temp)
    for k = 1 : 3
        tempROIsNb=find([temp{k}]<=Numbers(idx+1));
        if tempROIsNb            
            ROIsNb=[ROIsNb ; temp{1,k}(tempROIsNb)];
            temp{1,k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb ; repmat(k,length(tempROIsNb),1)];
        end
    end
    if ROIsNb
        imagename=regexp(filename,'_output_analysis_matlab.mat','split');
        imagename=strcat('AVG_',imagename{1},'.tif');
        image=double(imread(strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\',imagename)));image=image/max(max(image));image=image*128;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        load(filename, 'ROIs');
        Raster=Spikes(ROIsNb,:);
        ROIsNb=ROIsNb-Numbers(idx);
        ROIs=ROIs(:,ROIsNb);
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=reshape(ROIs(:,k),size(image,1),size(image,2));            
            image2=image2+ROI;image2=(image2/max(max(image2)))*200;image2=uint8(image2);
            image3(:,:,ClusterNb(k))=image3(:,:,ClusterNb(k))+image2;
        end
        %image3(:,:,3)=image;
    end
    %image3=uint8(image3);
    name=strcat('Spike_Kmeans',imagename(4:end));
    imwrite(image3,name,'tif');
    figure('Visible','Off');subplot(1,2,1);subimage(image3);subplot(1,2,2);imagesc(Raster,[0 5]);colormap jet;
    print(strcat(name(1:end-4),'.emf'),'-dmeta');
end

clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster

figure;hold on;colorstring='rgb'
counter=1;
for i=Select_model_Clusters
    subplot(3,1,counter);plot(Cmap_selectZS(i,:)','Color',colorstring(counter));
    counter=counter+1;
end
hold off;
print('Model_KMeansPlot.emf','-dmeta');

if Test
plot(Test);
end

Selected_idx=temp{1,1};
Selected_idx=[Selected_idx ; temp{1,2}];
Selected_idx=[Selected_idx ; temp{1,3}];
Selected=Calcium(Selected_idx,:);
Selected_sp=Spikes(Selected_idx,:);

figure;
for i=1:3    
    subplot(3,1,i);plot(select_Cmap(i,:),'Color',colors(i,:));
end

print('KMeansPlot_GoodBetas.emf','-dmeta');

GoodClustersData=[];
for i=1:length(Goodbetas_Select)
    GoodClustersData(i).FluoTraces=SelectCorr_Thr5(idxKmeans_DF==GoodBetas(i),:);
    GoodClustersData(i).DF=Select_DF_Thr5(idxKmeans_DF==GoodBetas(i),:)*100;
    GoodClustersData(i).Mean=mean(GoodClustersData(i).DF,1);
    GoodClustersData(i).STD=std(GoodClustersData(i).DF,1,1);    
end

Mean_GB=[GoodClustersData.Mean];Mean_GB=reshape(Mean_GB,length(GoodClustersData(1).Mean),length(GoodClustersData));
STD_GB=[GoodClustersData.STD];STD_GB=reshape(STD_GB,length(GoodClustersData(1).STD),length(GoodClustersData));

nb_GB=ones(length(GoodClustersData),length(GoodClustersData(1).Mean));
for i=1:length(GoodBetas)
    nb_GB(i,:)=size(GoodClustersData(i).DF,1)*nb_GB(i,:);
end
nb_GB=nb_GB';

files_info=[];files_info.names=sorted_files;files_info.index=sorted_index2;
ClustersData=[];
counter=1;counter2=1;
for i=GoodBetas
    temp=find(idx_final==i);
    for i=1:length(files_info.index)
        idx_clust=find(temp<files_info.index(i));
        temp(idx_clust)=[];
        [Fish,~]=regexp(files_info.names(i),'F(\d+)','tokens','match');
        [Plane,~]=regexp(files_info.names(i),'(\d+)um','tokens','match');
        if length(idx_clust)
            ClustersData{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=length(idx_clust);
        else
            ClustersData{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=0;
        end
    end
    counter=counter+1;
end
idxempty=cellfun('isempty',ClustersData);
ClustersData(idxempty)={0};
ClustersData=cell2mat(ClustersData);
clearvars counter i temp idx_clust Fish Plane idxempty counter2
FishData=squeeze(sum(ClustersData,2));
PlaneData=squeeze(sum(ClustersData,1));

OverallData=[];
counter=1;
for i=1:length(files_info.index)		
	[Fish,~]=regexp(files_info.names(i),'F(\d)','tokens','match');
	[Plane,~]=regexp(files_info.names(i),'(\d+)um','tokens','match');
    if i==1
        OverallData{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=files_info.index(i);    
    else
        OverallData{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=files_info.index(i)-files_info.index(i-1);
    end
end
idxempty=cellfun('isempty',OverallData);
OverallData(idxempty)={0};
OverallData=cell2mat(OverallData);
clearvars counter i temp idx_clust Fish Plane idxempty

FishOverallData=squeeze(sum(OverallData,2));
PlaneOverallData=squeeze(sum(OverallData,1));

NormalizedFishData = bsxfun(@rdivide,FishData,sum(FishData));
NormalizedPlaneData = bsxfun(@rdivide,PlaneData,sum(PlaneData));

Goodbetas_Select=Goodbetas([1,4,7]);
GoodClustersData_Select=GoodClustersData([1 4 7]);


DF_orderedbyClusters=[];
DF_orderedbyClusters=DF(idxKmeans_DF==1,:);
for i=2:max(idxKmeans_DF)
    temp=DF(idxKmeans_DF==i,:);
    DF_orderedbyClusters=vertcat(DF_orderedbyClusters,temp);
end
clearvars temp

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1049, 895]);
imagesc(DF_orderedbyClusters,[0 20]);colormap(jet);

Random2k_DF=datasample(DF_orderedbyClusters,2000);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1049, 895]);
imagesc(Random2k_DF,[0 20]);colormap(jet);

DF_orderedbyCorr=[];
DF_orderedbyCorr=DF(Nb_corr_selected==4,:);
for i=5:max(Nb_corr_selected)
    temp=DF(Nb_corr_selected==i,:);
    DF_orderedbyCorr=vertcat(DF_orderedbyCorr,temp);
end
clearvars temp

FishPlaneClusterData=[];
counter=1;
for j=Goodbetas_Select
    temp=find(idx_final==j);
    temp2=find(idxKmeans_DF==j);
    for i=1:length(files_info.index)
        idx_clust=find(temp<files_info.index(i));
        tempDF=DF(temp2(idx_clust),:);
        temp(idx_clust)=[];
        temp2(idx_clust)=[];
        [Fish,~]=regexp(files_info.names(i),'F(\d)','tokens','match');
        [Plane,~]=regexp(files_info.names(i),'(\d+)um','tokens','match');
        if isempty(idx_clust)
            FishPlaneClusterData{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=0;            
        else
            FishPlaneClusterData{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=tempDF;
        end
    end
    counter=counter+1;
end
idxempty=cellfun('isempty',FishPlaneClusterData);
FishPlaneClusterData(idxempty)={0};
clearvars counter i j k temp temp2 tempDF idx_clust Fish Plane idxempty

ClusterDF_Fish=[];
temp=[];
for i=1:3
    for j=1:size(FishPlaneClusterData,1)
        counter=1;
        for k=1:size(FishPlaneClusterData,2)
            if ~FishPlaneClusterData{j,k,i}==0
                if counter
                    temp=FishPlaneClusterData{j,k,i};counter=0;
                else
                    temp=vertcat(temp,FishPlaneClusterData{j,k,i});                    
                end
            end
        end
        ClusterDF_Fish{j,i}=temp;
    end    
end
clearvars counter i j k temp temp2 tempDF idx_clust Fish Plane idxempty

MeanClusterDF_Fish=zeros(size(ClusterDF_Fish,1),size(ClusterDF_Fish,2),size(DF,2));
for i=1:size(ClusterDF_Fish,1)
    for j=1:size(ClusterDF_Fish,2)
        MeanClusterDF_Fish(i,j,:)=mean(ClusterDF_Fish{i,j},1);
    end
end
clearvars  i j 

for i=1:size(MeanClusterDF_Fish,1)
    temp=squeeze(MeanClusterDF_Fish(i,:,:))';
    filename=strcat('Mean_Goodclust_fish_',num2str(i));
    save(filename,'temp');
end
clearvars i temp filename
    

ClustersData_Select=ClustersData(:,:,[1 4 7]);

counter=1;counter2=1;
idx_final=zeros(length(idx_corr),1);
for i=1:length(idx_corr)
    if idx_corr(i)
        if any(counter==idx_Threshold_5to200)    
            idx_final(i)=idxKmeans_DF(counter2);
            counter=counter+1;counter2=counter2+1;
        else
            counter=counter+1;
        end
    end
end
clearvars counter counter2;

DF=DeltaF2(double(SelectCorr),21,11);
max_DF=max(DF,[],2);
min_DF=min(DF,[],2);
idx_Threshold_5to200=find(max_DF>0.05 & max_DF<2 & min_DF>-0.1);
Select_DF_Thr5=DF(idx_Threshold_5to200,:);

random_5k=datasample(Select_DF_Thr5,5000);
eva = evalclusters(random_5k,'kmeans','silhouette','Distance','correlation','KList',[1:100]);
options = statset('UseParallel',1); [idxKmeans Cmap]=kmeans(Select_DF_Thr5,eva.OptimalK,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');

WholeBrainClusters_corr=zeros(size(Select_DF_Thr5,1),3);
parfor i=1:size(Select_DF_Thr5,1)
    for k=1:3
        temp=corrcoef(WholeBrainClusters(:,k),Select_DF_Thr5(i,:));
        temp=temp(1,2);
        WholeBrainClusters_corr(i,k)=temp;
    end    
end

[Max_corr_clusters(:,1), Max_corr_clusters(:,2)]=max(WholeBrainClusters_corr,[],2);idx_Max_corrClust=find(Max_corr_clusters(:,1)>0.5);
Select_corr_clusters=Max_corr_clusters(idx_Max_corrClust,:);

[sort_corr, idx_sort_corr]=sortrows(Nb_corr(idx_corr(idx_Threshold_5to200)));
Sort_Select_DF_Thr5=Select_DF_Thr5(idx_sort_corr,:);

STD_sp=std(Sp_infer,1,2);
histogram(STD_sp)
idx_STD_bool=STD_sp>20;idx_STD=find(STD_sp>20);
Select_Sp=Sp_infer(idx_STD,:);
random_5k=datasample(Select_Sp,5000);
eva = evalclusters(random_5k,'kmeans','silhouette','Distance','correlation','KList',[1:100]);
options = statset('UseParallel',1); [idxKmeans Cmap]=kmeans(Select_Sp,eva.OptimalK,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');

index = find(ismember(TransfoMatricesFishPlane, strcat(num2str(Fish),'_',num2str(Plane))));

temp=[TransfoMatricescleared{:,2}];

x_offset=640;y_offset=480;
temp_x=475;temp_y=419;
temp_x=temp_x-x_offset;temp_y=temp_y-y_offset;
angle=degtorad(8.8);
temp_xb=(temp_x*cos(angle)-temp_y*sin(angle))+x_offset;
temp_yb=(temp_y*cos(angle)+temp_x*sin(angle))+y_offset;

F5_centroids=Centroids_CTRL(49779:52468,:);
x_offset=640;y_offset=480;
figure;scatter(F5_centroids(:,1),F5_centroids(:,2));

angle=degtorad(-91);
F5_centroids_corrected(:,1)=F5_centroids(:,1)*cos(angle)-F5_centroids(:,2)*sin(angle);
F5_centroids_corrected(:,2)=F5_centroids(:,1)*sin(angle)+F5_centroids(:,2)*cos(angle);
figure;scatter(F5_centroids_corrected(:,1),F5_centroids_corrected(:,2));
F5_centroids_corrected(:,1)=F5_centroids_corrected(:,1)+abs(min(F5_centroids_corrected(:,1)));
F5_centroids_corrected(:,2)=F5_centroids_corrected(:,2)+abs(min(F5_centroids_corrected(:,2)));
figure;scatter(F5_centroids_corrected(:,1),F5_centroids_corrected(:,2));

index = find(ismember(TransfoMatricesFishPlane, '5_100'));Fish=5;
F5_centroids_corrected(:,1)=F5_centroids_corrected(:,1)+abs(TransfoMatricescleared{index,3}*cos(angles(Fish))-TransfoMatricescleared{index,2}*sin(angles(Fish)));
F5_centroids_corrected(:,2)=F5_centroids_corrected(:,2)+abs(TransfoMatricescleared{index,3}*cos(angles(Fish))+TransfoMatricescleared{index,2}*sin(angles(Fish)));

index = find(ismember(TransfoMatricesFishPlane, '1_100'));Fish=1;
F1_centroids_corrected(:,1)=F1_centroids(:,1)+abs(TransfoMatricescleared{index,2}*cos(angles(Fish))-TransfoMatricescleared{index,3}*sin(angles(Fish)));
F1_centroids_corrected(:,2)=F1_centroids(:,2)+abs(TransfoMatricescleared{index,2}*cos(angles(Fish))+TransfoMatricescleared{index,3}*sin(angles(Fish)));
figure;scatter(F1_centroids_corrected(:,1),F1_centroids_corrected(:,2));
