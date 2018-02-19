idx_corr_CTRL=find(Nb_corr_CTRL>3);
idx_corr_Neo=find(Nb_corr_Neo>3);
STD_CTRL=std(Sp_infer_CTRL,1,2);STD_Neo=std(Sp_infer_Neo,1,2);
ALL_Spikes=[Sp_infer_CTRL(STD_CTRL>2,:);Sp_infer_Neo(STD_Neo>2,:)];
options = statset('UseParallel',1); [idxKmeans Cmap]=kmeans(ALL_Spikes,100,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
[Model_DF,GoodBetas]=Test_Regress(Cmap,aud8freq_sp,idxKmeans,0.5);
BadBetas=find(rsq <0.7 & rsq >0.3);

ALL_DF=vertcat(DF,Select_DF_Thr5);
%DF == CTRL
%Select_DF_Thr5 == Neo
idx_final_CTRL=idx_corr(idx_Threshold_5to200);


GoodClustersData=[];
for i=1:length(GoodBetas)    
    GoodClustersData(i).DF=ALL_DF(idxKmeans==GoodBetas(i),:)*100;
    GoodClustersData(i).Mean=mean(GoodClustersData(i).DF,1);
    GoodClustersData(i).STD=std(GoodClustersData(i).DF,1,1);    
end

GoodClustersData_CTRL=[];
for i=1:length(GoodBetas)    
    GoodClustersData_CTRL(i).DF=ALL_DF(idxKmeans(1:length(DF))==GoodBetas(i),:)*100;
    GoodClustersData_CTRL(i).Mean=mean(GoodClustersData_CTRL(i).DF,1);
    GoodClustersData_CTRL(i).STD=std(GoodClustersData_CTRL(i).DF,1,1);    
end

BadClustersData=[];
for i=1:length(BadBetas)    
    BadClustersData(i).DF=ALL_DF(idxKmeans==BadBetas(i),:)*100;
    BadClustersData(i).Mean=mean(BadClustersData(i).DF,1);
    BadClustersData(i).STD=std(BadClustersData(i).DF,1,1);    
end

Mean_GB=[GoodClustersData.Mean];Mean_GB=reshape(Mean_GB,length(GoodClustersData(1).Mean),length(GoodClustersData));
STD_GB=[GoodClustersData.STD];STD_GB=reshape(STD_GB,length(GoodClustersData(1).STD),length(GoodClustersData));
nb_GB=ones(length(GoodClustersData),length(GoodClustersData(1).Mean));
for i=1:length(GoodBetas)
    nb_GB(i,:)=size(GoodClustersData(i).DF,1)*nb_GB(i,:);
end
nb_GB=nb_GB';

BadMean_GB=[BadClustersData.Mean];BadMean_GB=reshape(BadMean_GB,length(BadClustersData(1).Mean),length(BadClustersData));
BadSTD_GB=[BadClustersData.STD];BadSTD_GB=reshape(BadSTD_GB,length(BadClustersData(1).STD),length(BadClustersData));
Badnb_GB=ones(length(BadClustersData),length(BadClustersData(1).Mean));
for i=1:length(BadBetas)
    Badnb_GB(i,:)=size(BadClustersData(i).DF,1)*Badnb_GB(i,:);
end
Badnb_GB=Badnb_GB';
clearvars i ans counter2

NeoOrCTRL=[];
for i=1:length(GoodBetas)    
    NeoOrCTRL{i,1}=find(idxKmeans(1:length(DF))==GoodBetas(i));
    NeoOrCTRL{i,2}=find(idxKmeans(length(DF)+1:end)==GoodBetas(i));
end

counter=1;
idx_final_CTRL=zeros(length(idx_corr),1);
for i=1:length(idx_corr)
    if idx_corr(i)
            idx_final_CTRL(i)=idxKmeans(counter);
            counter=counter+1;        
    end
end
clearvars counter counter2;

counter=1;counter2=1;
idx_final_Neo=zeros(length(idx_corr_bool),1);
for i=1:length(idx_corr_bool)
    if idx_corr_bool(i)
        if ~idx_Threshold_5to200_bool(counter)
            counter=counter+1;
        else
            idx_final_Neo(i)=idxKmeans(length(DF)+counter2);
            counter=counter+1;counter2=counter2+1;
        end
    end
end
clearvars counter counter2;

ClustersData_CTRL=[];
counter=1;counter2=1;
for i=GoodBetas
    temp=find(idx_final_CTRL==i);    
    for i=1:length(files_info_CTRL.index)
        idx_clust=find(temp<files_info_CTRL.index(i));
        temp(idx_clust)=[];
        [Fish,~]=regexp(files_info_CTRL.names(i),'F(\d+)','tokens','match');
        [Plane,~]=regexp(files_info_CTRL.names(i),'(\d+)um','tokens','match');
        if length(idx_clust)
            ClustersData_CTRL{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=length(idx_clust);
        else
            ClustersData_CTRL{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=0;
        end
    end
    counter=counter+1;
end
idxempty=cellfun('isempty',ClustersData_CTRL);
ClustersData_CTRL(idxempty)={0};
ClustersData_CTRL=cell2mat(ClustersData_CTRL);
clearvars counter i temp idx_clust Fish Plane idxempty counter2
FishData_CTRL=squeeze(sum(ClustersData_CTRL,2));
PlaneData_CTRL=squeeze(sum(ClustersData_CTRL,1));

ClustersData_Neo=[];
counter=1;
for j=GoodBetas
    temp=find(idx_final_Neo==j);
    [Fish,~]=regexp(files_info_Neo.names(1),'F(\d+)','end');Fish=Fish{1,1};
    tempcomp=files_info_Neo.names(1);tempcomp=tempcomp{1,1}(1:Fish);
    counter2=1;
    for i=1:length(files_info_Neo.index)
        idx_clust=find(temp<files_info_Neo.index(i));
        temp(idx_clust)=[];
        [Fish,~]=regexp(files_info_Neo.names(i),'F(\d+)','end');Fish=Fish{1,1};
        tempcomp2=files_info_Neo.names(i);tempcomp2=tempcomp2{1,1}(1:Fish);        
        [Plane,~]=regexp(files_info_Neo.names(i),'(\d+)um','tokens','match');
        if ~strcmp(tempcomp,tempcomp2)
            tempcomp=tempcomp2;
            counter2=counter2+1;
        end
        if length(idx_clust)
            ClustersData_Neo{counter2,1+str2num(Plane{1,1}{1}{1})/10,counter}=length(idx_clust);
        else
            ClustersData_Neo{counter2,1+str2num(Plane{1,1}{1}{1})/10,counter}=0;
        end
    end
    counter=counter+1;
end
idxempty=cellfun('isempty',ClustersData_Neo);
ClustersData_Neo(idxempty)={0};
ClustersData_Neo=cell2mat(ClustersData_Neo);
clearvars counter i temp idx_clust Fish Plane idxempty
FishData_Neo=squeeze(sum(ClustersData_Neo,2));
PlaneData_Neo=squeeze(sum(ClustersData_Neo,1));

OverallData_CTRL=[];
counter=1;
for i=1:length(files_info_CTRL.index)		
	[Fish,~]=regexp(files_info_CTRL.names(i),'F(\d)','tokens','match');
	[Plane,~]=regexp(files_info_CTRL.names(i),'(\d+)um','tokens','match');
    if i==1
        OverallData_CTRL{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=files_info_CTRL.index(i);    
    else
        OverallData_CTRL{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=files_info_CTRL.index(i)-files_info_CTRL.index(i-1);
    end
end
idxempty=cellfun('isempty',OverallData_CTRL);
OverallData_CTRL(idxempty)={0};
OverallData_CTRL=cell2mat(OverallData_CTRL);
clearvars counter i temp idx_clust Fish Plane idxempty

FishOverallData_CTRL=squeeze(sum(OverallData_CTRL,2));
PlaneOverallData_CTRL=squeeze(sum(OverallData_CTRL,1));

OverallData_Neo=[];
counter=1;counter2=1;tempcomp=files_info_Neo.names(1);[Fish,~]=regexp(files_info_Neo.names(1),'F(\d+)','end');Fish=Fish{1,1};tempcomp=tempcomp{1,1}(1:Fish);
for i=1:length(files_info_Neo.index)		
	[Fish,~]=regexp(files_info_Neo.names(i),'F(\d+)','end');Fish=Fish{1,1};
    tempcomp2=files_info_Neo.names(i);tempcomp2=tempcomp2{1,1}(1:Fish);        
	[Plane,~]=regexp(files_info_Neo.names(i),'(\d+)um','tokens','match');
	if ~strcmp(tempcomp,tempcomp2)
        tempcomp=tempcomp2;
        counter2=counter2+1;
    end
	
    if i==1
        OverallData_Neo{counter2,1+str2num(Plane{1}{1}{1})/10,counter}=files_info_Neo.index(i);    
    else
        OverallData_Neo{counter2,1+str2num(Plane{1}{1}{1})/10,counter}=files_info_Neo.index(i)-files_info_Neo.index(i-1);
    end
end
idxempty=cellfun('isempty',OverallData_Neo);
OverallData_Neo(idxempty)={0};
OverallData_Neo=cell2mat(OverallData_Neo);
clearvars counter i temp idx_clust Fish Plane idxempty tempcomp tempcomp2 test

FishOverallData_Neo=squeeze(sum(OverallData_Neo,2));
PlaneOverallData_Neo=squeeze(sum(OverallData_Neo,1));

NormalizedFishData = bsxfun(@rdivide,vertcat(FishData_CTRL,FishData_Neo),sum(vertcat(FishData_CTRL,FishData_Neo)));
NormalizedPlaneData = bsxfun(@rdivide,vertcat(PlaneData_CTRL,PlaneData_Neo),sum(vertcat(PlaneData_CTRL,PlaneData_Neo)));

PlaneOverallData=PlaneOverallData_CTRL(1:23)+PlaneOverallData_Neo;
PlaneOverallData=[PlaneOverallData PlaneOverallData_CTRL(24:end)];
NormalizedPlaneOverallData=bsxfun(@rdivide,PlaneOverallData',sum(PlaneOverallData'));

temp=[Mean_GB(:,1) STD_GB(:,1) nb_GB(:,1)];
for i=2:length(GoodBetas)
    temp=horzcat(temp,[Mean_GB(:,i) STD_GB(:,i) nb_GB(:,i)]);
end

temp=[BadMean_GB(:,1) BadSTD_GB(:,1) Badnb_GB(:,1)];
for i=2:length(BadBetas)
    temp=horzcat(temp,[BadMean_GB(:,i) BadSTD_GB(:,i) Badnb_GB(:,i)]);
end

CTRL_Mask=[];
for i=1:length(files_info_CTRL.names)
    filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f_final\Results\','Maskb_',files_info_CTRL.names(i));filename=filename{1,1};
    CTRL_Mask{i}=imread(filename);    
end
clearvars i idxsort_files filename;

Neo_Mask=[];
for i=1:length(files_info_Neo.names)
    filename=strcat('D:\Pictures\processed\Tonotropy\neomycin\','Maskb_',files_info_Neo.names(i));filename=filename{1,1};
    Neo_Mask{i}=imread(filename);    
end
clearvars i idxsort_files filename;

colors = distinguishable_colors(length(GoodBetasSelect),[1 1 1 ; 0 0 0]);
colors = colors*256;
for i=1:length(CTRL_Mask)
    B=CTRL_Mask{i};B=double(B);B(B==0)=NaN;
    filename=files_info_CTRL.names(i);
    filename=strcat('D:\Pictures\processed\Tonotropy\AVG_Cropped_CTRL\AVG_',filename{1});
    image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(max((image4)));image4=image4*128;image4=double(repmat(image4,1,1,3));image4=uint8(image4);
    %image=zeros(size(B,1),size(B,2),3);
    for j=min(min(B)):max(max(B))
        if any(idx_final_CTRL(j)==GoodBetasSelect)
            clust=find(idx_final_CTRL(j)==GoodBetasSelect);
            for k=1:3              
              image2=image4(:,:,k);
              image2(B==j)=colors(clust,k);
              image4(:,:,k)=image2;
            end
        end
    end
    %image=uint8(image);
    filename=files_info_CTRL.names(i);
    name=strcat('Neo_CTRL_',filename{1});
    %filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\AVG_',sorted_files{i});
    %image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(prctile(image4,80));image4=image4*200;
    %image4=double(repmat(image4,1,1,3));image4=uint8(image4)+image;
    image4=uint8(image4);
    imwrite(image4,name,'tif');
end
clearvars image image2 image4 i j k name filename B clust;

colors = distinguishable_colors(length(GoodBetasSelect),[1 1 1 ; 0 0 0]);
colors = colors*256;
for i=1:length(Neo_Mask)
    B=Neo_Mask{i};B=double(B);B(B==0)=NaN;
    filename=files_info_Neo.names(i);
    filename=strcat('D:\Pictures\processed\Tonotropy\neomycin\AVG_',filename{1});
    image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(max((image4)));image4=image4*128;image4=double(repmat(image4,1,1,3));image4=uint8(image4);
    %image=zeros(size(B,1),size(B,2),3);
    for j=min(min(B)):max(max(B))
        if any(idx_final_Neo(j)==GoodBetasSelect)
            clust=find(idx_final_Neo(j)==GoodBetasSelect);
            for k=1:3              
              image2=image4(:,:,k);
              image2(B==j)=colors(clust,k);
              image4(:,:,k)=image2;
            end
        end
    end
    %image=uint8(image);
    filename=files_info_Neo.names(i);
    name=strcat('Neo_CTRL_',filename{1});
    %filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\AVG_',sorted_files{i});
    %image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(prctile(image4,80));image4=image4*200;
    %image4=double(repmat(image4,1,1,3));image4=uint8(image4)+image;
    image4=uint8(image4);
    imwrite(image4,name,'tif');
end
clearvars image image2 image4 i j k name filename B clust;
clearvars ans angle i j k counter counter2 Fish Plane

BadClustersData_CTRL=[];
counter=1;counter2=1;
for i=BadBetas
    temp=find(idx_final_CTRL==i);    
    for i=1:length(files_info_CTRL.index)
        idx_clust=find(temp<files_info_CTRL.index(i));
        temp(idx_clust)=[];
        [Fish,~]=regexp(files_info_CTRL.names(i),'F(\d+)','tokens','match');
        [Plane,~]=regexp(files_info_CTRL.names(i),'(\d+)um','tokens','match');
        if length(idx_clust)
            BadClustersData_CTRL{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=length(idx_clust);
        else
            BadClustersData_CTRL{str2num(Fish{1}{1}{1}),1+str2num(Plane{1}{1}{1})/10,counter}=0;
        end
    end
    counter=counter+1;
end
idxempty=cellfun('isempty',BadClustersData_CTRL);
BadClustersData_CTRL(idxempty)={0};
BadClustersData_CTRL=cell2mat(BadClustersData_CTRL);
clearvars counter i temp idx_clust Fish Plane idxempty counter2
BadFishData_CTRL=squeeze(sum(BadClustersData_CTRL,2));
BadPlaneData_CTRL=squeeze(sum(BadClustersData_CTRL,1));

BadClustersData_Neo=[];
counter=1;
for j=BadBetas
    temp=find(idx_final_Neo==j);
    [Fish,~]=regexp(files_info_Neo.names(1),'F(\d+)','end');Fish=Fish{1,1};
    tempcomp=files_info_Neo.names(1);tempcomp=tempcomp{1,1}(1:Fish);
    counter2=1;
    for i=1:length(files_info_Neo.index)
        idx_clust=find(temp<files_info_Neo.index(i));
        temp(idx_clust)=[];
        [Fish,~]=regexp(files_info_Neo.names(i),'F(\d+)','end');Fish=Fish{1,1};
        tempcomp2=files_info_Neo.names(i);tempcomp2=tempcomp2{1,1}(1:Fish);        
        [Plane,~]=regexp(files_info_Neo.names(i),'(\d+)um','tokens','match');
        if ~strcmp(tempcomp,tempcomp2)
            tempcomp=tempcomp2;
            counter2=counter2+1;
        end
        if length(idx_clust)
            BadClustersData_Neo{counter2,1+str2num(Plane{1,1}{1}{1})/10,counter}=length(idx_clust);
        else
            BadClustersData_Neo{counter2,1+str2num(Plane{1,1}{1}{1})/10,counter}=0;
        end
    end
    counter=counter+1;
end
idxempty=cellfun('isempty',BadClustersData_Neo);
BadClustersData_Neo(idxempty)={0};
BadClustersData_Neo=cell2mat(BadClustersData_Neo);
clearvars counter i temp idx_clust Fish Plane idxempty
BadFishData_Neo=squeeze(sum(BadClustersData_Neo,2));
BadPlaneData_Neo=squeeze(sum(BadClustersData_Neo,1));

BadNormalizedFishData = bsxfun(@rdivide,vertcat(BadFishData_CTRL,BadFishData_Neo),sum(vertcat(BadFishData_CTRL,BadFishData_Neo)));
BadNormalizedPlaneData = bsxfun(@rdivide,vertcat(BadPlaneData_CTRL,BadPlaneData_Neo),sum(vertcat(BadPlaneData_CTRL,BadPlaneData_Neo)));

figure;
for i=1:4
    subplot(4,1,i);plot(Cmap(GoodBetasSelect(i),:),'Color',colors(i,:)/256);
end

Nb_corr_Neo=load('neomycin-FluoTraces-Maskb-SelectCorr.mat', 'Nb_corr')
Nb_corr_Neo(idx_corr_bool);ans(idx_Threshold_5to200_bool);
Nb_corr_Neo=ans;
Nb_corr_CTRL=load('GCaMP6f2-FluoTraces-Maskb-SelectCorr.mat', 'Nb_corr')
Nb_corr_CTRL(idx_corr);Nb_corr_CTRL=ans;
Nb_corr=vertcat(Nb_corr_CTRL,Nb_corr_Neo);
Nb_corr=vertcat(Nb_corr_CTRL,Nb_corr_Neo);
[sort_corr, idx_sort_corr]=sortrows(Nb_corr);
ALL_DF=ALL_DF(idx_sort_corr,:);

x = linspace(0.2,140,700);
y = linspace(1,size(ALL_DF,1),450000);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ALL_DF*100,[0 20]);colormap(jet);

y = linspace(1,size(Cmap,1),size(Cmap,1));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,Cmap*100,[0 20]);colormap(jet);

parfor i=1:size(ALL_DF,1)
    mdl=stepwiselm(aud8freq',ALL_DF(i,:),'linear','Criterion','adjrsquared','Upper','interactions','verbose',0);
    model_DF_Thr5(i).coef=mdl.Coefficients;
    model_DF_Thr5(i).MSE=mdl.MSE;
    model_DF_Thr5(i).Fitted=mdl.Fitted;
    model_DF_Thr5(i).rsquared=mdl.Rsquared.Adjusted;
end

coefficients_above02=[];
idx_02=find(rsq>0.2);
for idx=1:length(idx_02)
    coef=[model_DF_Thr5(idx_02(idx)).coef];
    for coef_idx=2:9
        if coef.pValue(coef_idx)<0.05
            coefficients_above02{idx,coef_idx-1}=coef.Estimate(coef_idx);
        end
    end    
end
idxempty=cellfun('isempty',coefficients_above02);
coefficients_above02(idxempty)={0};
clearvars idxempty idx coef_idx coef
coefficients_above02=cell2mat(coefficients_above02);
[max_coef,idx_max]=max(coefficients_above02,[],2);

%import TransfoMatrices as cell
uiopen('D:\Pictures\processed\Tonotropy\GCaMP6f_final\NeoandCTRL_reconstruct\TransfoMatrices_cleared.txt',1)
%1st row is Fish_Plane, 2nd row is x translation, 3rd row is y translation
TransfoMatricesFishPlane=TransfoMatricescleared(:,1);

CenteringMatrix=zeros(8,26,2);
list=regexp(Centereverything(:,1),'_','split');
for i=1:length(list)
    Fish=str2num(list{i}{1});Plane=(str2num(list{i}{2})-40)/10;    
    CenteringMatrix(Fish,Plane,1)=Centereverything{i,2};CenteringMatrix(Fish,Plane,2)=Centereverything{i,3};
end
clearvars i list Fish Plane j counter2 ans index

for Fish=1:size(CenteringMatrix,1)
    %if Fish~=5
    indices=find(CenteringMatrix(Fish,:,1)~=0 | CenteringMatrix(Fish,:,2)~=0);
    counter=2;counter2=1;
    for Plane=2:size(CenteringMatrix,2)
        if Plane==indices(counter)
            if length(indices)>=counter+1
                counter2=counter;
                counter=counter+1;
            else
                CenteringMatrix(Fish,indices(counter)+1:end,1)=CenteringMatrix(Fish,indices(counter),1);
                CenteringMatrix(Fish,indices(counter)+1:end,2)=CenteringMatrix(Fish,indices(counter),2);
                break
            end
        else
            CenteringMatrix(Fish,Plane,1)=(((indices(counter)-Plane)*CenteringMatrix(Fish,indices(counter2),1))+((Plane-indices(counter2))*CenteringMatrix(Fish,indices(counter),1)))/(indices(counter)-indices(counter2));
            CenteringMatrix(Fish,Plane,2)=(((indices(counter)-Plane)*CenteringMatrix(Fish,indices(counter2),2))+((Plane-indices(counter2))*CenteringMatrix(Fish,indices(counter),2)))/(indices(counter)-indices(counter2));
            %CenteringMatrix(Fish,Plane,1)=CenteringMatrix(Fish,indices(counter2),1)+CenteringMatrix(Fish,indices(counter),1);
        end
    end
    %end
end

Plane_correction=[-10 50 0 -10 20 10 -20 0];
%Plane_correction=[0 0 0 0 0 0 0 0];

%NEW VERSION

angles=degtorad([-4.4 -8.3 -1.5 -1.6 88 98 -2 -2.5]);
Centroids_CTRL_corrected=zeros(size(Centroids_CTRL,1),4);
%for i=1:length(files_info_CTRL.index)
%for i=[1:1:93 141:1:190]
for i=[1:1:72 141:1:190]
%for i=[141:1:190]
    if i==1
        counter=1;
    else
        counter=files_info_CTRL.index(i-1)+1;
    end
	[Fish,~]=regexp(files_info_CTRL.names(i),'F(\d)','tokens','match');Fish=str2num(Fish{1}{1}{1});
    [Plane,~]=regexp(files_info_CTRL.names(i),'(\d+)um','tokens','match');Plane=str2num(Plane{1}{1}{1});
    idx=files_info_CTRL.index(i);  
    temp_x=Centroids_CTRL(counter:idx,2);
    temp_y=Centroids_CTRL(counter:idx,1);
    index = find(ismember(TransfoMatricesFishPlane, strcat(num2str(Fish),'_',num2str(Plane))));index_plane=(Plane-40)/10;
    if index_plane<1
    	index_plane=1;
    end    
    temp_x = temp_x+TransfoMatricescleared{index,2}+CenteringMatrix(Fish,index_plane,1);
    temp_y = temp_y+TransfoMatricescleared{index,3}+CenteringMatrix(Fish,index_plane,2);
	temp_xb=temp_x*cos(angles(Fish))-temp_y*sin(angles(Fish));%temp_xb=temp_xb+abs(min(temp_xb));
	temp_yb=temp_y*cos(angles(Fish))+temp_x*sin(angles(Fish));%temp_yb=temp_yb+abs(min(temp_yb));
    if min(temp_xb)<0
        temp_xb=temp_xb+1080;
    end
    if min(temp_yb)<0
        temp_yb=temp_yb+1280;
    end
    Centroids_CTRL_corrected(counter:idx,1)=temp_xb;
    Centroids_CTRL_corrected(counter:idx,2)=temp_yb;        
    Centroids_CTRL_corrected(counter:idx,3)=-Plane+Plane_correction(Fish);
    Centroids_CTRL_corrected(counter:idx,4)=Fish;
    counter=idx+1;
end
%clearvars idx counter temp_x temp_y Fish Plane i

Centroids_clust_CTRL=[];
Centroids_index_CTRL=[];
Traces_clust_CTRL=[];
for i=1:4
    %Centroids_clust_CTRL{i}=Centroids_CTRL(find(idx_final_CTRL==GoodBetasSelect(i)),:);
    Centroids_clust_CTRL{i}=Centroids_CTRL_corrected(find(idx_final_CTRL==GoodBetasSelect(i)),:);
    if i==1
        Traces_clust_CTRL{i}=GoodClustersData_CTRL(i).DF;
    else
        Traces_clust_CTRL{i}=GoodClustersData_CTRL(i+1).DF;
    end
    Centroids_index_CTRL{i}=find(idx_final_CTRL==GoodBetasSelect(i));
end
clearvars i

temp=vertcat(Centroids_clust_CTRL{:});
temp2=zeros(size(temp,1),1);
temp3=vertcat(Traces_clust_CTRL{:});
counter=1;
for i=1:4
    temp2(counter:counter+size(Centroids_clust_CTRL{i}-1,1))=i;
    counter=counter+size(Centroids_clust_CTRL{i},1);
end
temp2(end)=[];
temp=[temp temp2];
idx=find(temp(:,1)==0);
temp(idx,:)=[];
temp3(idx,:)=[];
Select_ClusterTraces=temp3;
Select_Clustercoord=temp;clearvars i temp2 counter temp temp3

figure;
for i=1:4
    subplot(2,2,i);scatter(Centroids_clust_CTRL{i}(:,1),Centroids_clust_CTRL{i}(:,2),6,colors(i,:)/256);
end



Images=zeros(1080,1280,length(unique(Select_Clustercoord(:,3))));
disk=strel('disk',3,0);disk=double(disk.getnhood);
correct=1+max(Select_Clustercoord(:,3))/10;
for j=1:4
    Images=zeros(1080,1280,length(unique(Select_Clustercoord(:,3))));
    test=find(Select_Clustercoord(:,5)==j);
	for i=1:length(test)
        Plane=(-Select_Clustercoord(test(i),3)/10)+correct;
        x=Select_Clustercoord(test(i),1);y=Select_Clustercoord(test(i),2);
        Images(x-3:x+3,y-3:y+3,Plane)=disk*200;
    end
    ImageClustCoor{j}=Images;
end

for j=1:4
    filename=strcat(num2str(j),'_Clusters_coord.tif');
    Images=ImageClustCoor{j};
    saveastiff(uint8(Images), filename);
end

ImageClustCoor{1};
saveastiff(uint8(ans), 'ClusterBlue.tif');
ImageClustCoor{2};
saveastiff(uint8(ans), 'ClusterRed.tif');
ImageClustCoor{3};
saveastiff(uint8(ans), 'ClusterGreen.tif');
ImageClustCoor{4};
saveastiff(uint8(ans), 'ClusterPink.tif');

clearvars i j Images correct test x Plane filename


%temp is concatenated cluster_coord of selected fish
%Select_Clustercoord=temp;
ganglia=find(Select_Clustercoord(:,1)>640 & Select_Clustercoord(:,2)>620 & Select_Clustercoord(:,2)<950);
ganglia_coord=Select_Clustercoord(ganglia,:);
ganglia_Traces=Select_ClusterTraces(ganglia,:);
ganglia_ClustInfo=[];
for i=1:4
    temp=find(ganglia_coord(:,5)==i);
    ganglia_ClustInfo(i).DF=ganglia_Traces(temp,:);
    ganglia_ClustInfo(i).coord=ganglia_coord(temp,:);
    ganglia_ClustInfo(i).mean=mean(ganglia_ClustInfo(i).DF,1);
    ganglia_ClustInfo(i).STD=std(ganglia_ClustInfo(i).DF,1,1);
end
clearvars ganglia ganglia_coord ganglia_Traces

figure;
for i=1:4
    temp=find(ganglia_ClustInfo(i).coord(:,5)==i);temp=ganglia_ClustInfo(i).coord(temp,:);
    subplot(2,2,i);scatter(temp(:,1),temp(:,2),6,colors(i,:)/256);
end

angles_ganglia=degtorad([-52 -52 -55 -49 -45 -46 -41 -40]);
ganglia_coord=vertcat(ganglia_ClustInfo.coord);
ganglia_coord_corrected=[ganglia_coord(:,1:3) ganglia_coord(:,5)];
ganglia_coord_corrected(:,1)=ganglia_coord_corrected(:,1)-640;ganglia_coord_corrected(:,2)=ganglia_coord_corrected(:,2)-620;
for i=1:8
    temp=find(ganglia_coord(:,4)==i);
    if temp
        temp_y=ganglia_coord_corrected(temp,2);temp_z=ganglia_coord_corrected(temp,3);
        ganglia_coord_corrected(temp,2)=cos(angles_ganglia(i))*temp_y+sin(angles_ganglia(i))*temp_z;
        ganglia_coord_corrected(temp,3)=cos(angles_ganglia(i))*temp_z-sin(angles_ganglia(i))*temp_y;
    end
end
clearvars ganglia_coord

angles_test=([-52 -52 -55 -49 -45 -46 -41 -40]);
Test_coord_corrected=[Select_Clustercoord(:,1:3) Select_Clustercoord(:,5)];
for i=1:8
    temp=find(Select_Clustercoord(:,4)==i);
    if temp
        temp_y=Test_coord_corrected(temp,2);temp_z=Test_coord_corrected(temp,3);
        Test_coord_corrected(temp,2)=cos(angles_test(i))*temp_y+sin(angles_test(i))*temp_z;
        Test_coord_corrected(temp,3)=cos(angles_test(i))*temp_z-sin(angles_test(i))*temp_y;
    end
end
clearvars ganglia_coord

Torus=find(Select_Clustercoord(:,1)>750 & Select_Clustercoord(:,2)>450 & Select_Clustercoord(:,2)<620);
Torus_coord=Select_Clustercoord(Torus,:);
Torus_Traces=Select_ClusterTraces(Torus,:);
for i=1:4
    temp=find(Torus_coord(:,5)==i);
    Torus_ClustInfo(i).DF=Torus_Traces(temp,:);
    Torus_ClustInfo(i).coord=Torus_coord(temp,:);
    Torus_ClustInfo(i).mean=mean(Torus_ClustInfo(i).DF,1);
    Torus_ClustInfo(i).STD=std(Torus_ClustInfo(i).DF,1,1);
end
clearvars ganglia Torus_coord Torus_Traces

figure;
for i=1:4
    temp=find(Torus_coord(:,5)==i);temp=Torus_coord(temp,:);
    subplot(2,2,i);scatter(temp(:,1),temp(:,2),6,colors(i,:)/256);
end
Torus_coord_corrected=[Torus_coord(:,1:3) Torus_coord(:,5)];
Torus_coord_corrected(:,1)=Torus_coord_corrected(:,1)-750;Torus_coord_corrected(:,2)=Torus_coord_corrected(:,2)-450;
Torus_coord_Fish=zeros(4,8,3);
for i=1:8
    temp=find(Torus_coord(:,4)==i);
    if temp
        temp=[Torus_coord_corrected(temp,:) Torus_coord(temp,5)];
        for j=1:4    
            temp2=find(temp(:,5)==j);
            Torus_coord_Fish(j,i,1)=mean(temp(temp2,1));
            Torus_coord_Fish(j,i,2)=mean(temp(temp2,2));
            Torus_coord_Fish(j,i,3)=mean(temp(temp2,3));
        end
    end
end
clearvars temp counter i j

ganglia_coord_Fish=zeros(4,8,3);
ganglia_coord_perFish=[];
ganglia_coord=vertcat(ganglia_ClustInfo.coord);
Z_coord_fish=[];
for i=1:8
    temp=find(ganglia_coord(:,4)==i);
    if temp
        temp=[ganglia_coord_corrected(temp,:) ganglia_coord(temp,5)];
        for j=1:4    
            temp2=find(temp(:,5)==j);
            if length(temp2)>10
                ganglia_coord_Fish(j,i,1)=mean(temp(temp2,1));
                ganglia_coord_perFish{i,j}=temp(temp2,1:3);
                Z_coord_fish{i,j}=temp(temp2,3);
                ganglia_coord_Fish(j,i,2)=mean(temp(temp2,2));
                ganglia_coord_Fish(j,i,3)=mean(temp(temp2,3));
            end
            
        end
    end
end
clearvars temp counter i j ganglia_coord

figure;
Z_distrib=zeros(size(Z_coord_fish,1),4,30);
for i=1:size(Z_coord_fish,1)
    for j=1:4
        if Z_coord_fish{i,j}
            h=histogram(Z_coord_fish{i,j},30);
            Z_distrib(i,j,:)=h.Values;
        end
        
    end
end

temp=zeros(700,12);
for i=1:4
    temp(:,1+(i-1)*3)=ganglia_ClustInfo(i).mean;
    temp(:,2+(i-1)*3)=ganglia_ClustInfo(i).STD;
    temp(:,3+(i-1)*3)=size(ganglia_ClustInfo(i).DF,1);
end

temp=zeros(700,12);
for i=1:4
    temp(:,1+(i-1)*3)=Torus_ClustInfo(i).mean;
    temp(:,2+(i-1)*3)=Torus_ClustInfo(i).STD;
    temp(:,3+(i-1)*3)=size(Torus_ClustInfo(i).DF,1);
end



%idx_final for Linreg Data
counter=0;
idx_final_CTRL=zeros(length(idx_corr),1);
for i=1:length(idx_corr)
    if idx_corr(i)    	
    	counter=counter+1;
        counter2=find(idx_05==counter);
        if counter2
            idx_final_CTRL(i)=idx_max(counter2);
        end
    end
end
clearvars counter counter2;

colors = distinguishable_colors(8,[1 1 1 ; 0 0 0]);
colors = colors*256;
for i=1:length(CTRL_Mask)
    B=CTRL_Mask{i};B=double(B);B(B==0)=NaN;
    filename=files_info_CTRL.names(i);
    filename=strcat('D:\Pictures\processed\Tonotropy\AVG_Cropped_CTRL\AVG_',filename{1});
    image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(max((image4)));image4=image4*128;image4=double(repmat(image4,1,1,3));image4=uint8(image4);
    %image=zeros(size(B,1),size(B,2),3);
    for j=min(min(B)):max(max(B))
        if idx_final_CTRL(j)            
            for k=1:3              
              image2=image4(:,:,k);
              image2(B==j)=colors(idx_final_CTRL(j),k);
              image4(:,:,k)=image2;
            end
        end
    end
    %image=uint8(image);
    filename=files_info_CTRL.names(i);
    name=strcat('NCcoef_',filename{1});
    %filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\AVG_',sorted_files{i});
    %image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(prctile(image4,80));image4=image4*200;
    %image4=double(repmat(image4,1,1,3));image4=uint8(image4)+image;
    image4=uint8(image4);
    imwrite(image4,name,'tif');
end
clearvars image image2 image4 i j k name filename B clust ans angle counter counter2 Fish idx Plane;
clearvars index index_plane indices midline midline_F5_6

coefficients_above05=[];
idx_05=find(rsq>0.5);
for idx=1:length(idx_05)
    coef=[model_DF_Thr5(idx_05(idx)).coef];
    for coef_idx=2:9
        if coef.pValue(coef_idx)<0.05
            coefficients_above05{idx,coef_idx-1}=coef.Estimate(coef_idx);
        end
    end    
end
idxempty=cellfun('isempty',coefficients_above05);
coefficients_above05(idxempty)={0};
clearvars idxempty idx coef_idx coef
coefficients_above05=cell2mat(coefficients_above05);
[max_coef_05,idx_max_05]=max(coefficients_above05,[],2);

OverallData=zeros(15,31);
OverallData(1:8,:)=OverallData_CTRL;OverallData(9:15,1:23)=OverallData_Neo;
NormalizedOverallData=bsxfun(@rdivide,OverallData',sum(OverallData'));
NormalizedOverallData(NormalizedOverallData==0)=NaN;

OverallDataNorm.mean=nanmean(NormalizedOverallData*100,2);
OverallDataNorm.std=nanstd(NormalizedOverallData*100,1,2);

PlaneData=zeros(15,31,4);
PlaneData(1:8,:,:)=ClustersData_CTRL(:,:,[1 3 4 5]);PlaneData(9:15,1:23,:)=ClustersData_Neo(:,:,[1 3 4 5]);
NormalizedPlaneData=zeros(31,15,4);
temp=squeeze(sum(PlaneData~=0,1));
tempb=PlaneData;
for j=1:31
    for i=1:4
        if temp(j,i)<3
           tempb(:,j,i)=zeros(15,1);
        end
    end
end
for i=1:4
    NormalizedPlaneData(:,:,i)=bsxfun(@rdivide,tempb(:,:,i)',sum(tempb(:,:,i)'));
end
NormalizedPlaneData(NormalizedPlaneData==0)=NaN;
PlaneDataNorm.mean=squeeze(nanmean(NormalizedPlaneData*100,2));
PlaneDataNorm.std=squeeze(nanstd(NormalizedPlaneData*100,1,2));
clearvars temp tempb

progressbar
ToMerge{1}=[];
counter=1;
Merged_DF=ALL_DF_02;
for i=1:size(ALL_DF_02,1)
    max_corr=1;
    while max_corr>0.75
        corr_temp=zeros(size(Merged_DF,1),1);
        parfor j=1:size(Merged_DF,1)
            if i~=j
                temp=corrcoef(Merged_DF(i,:), Merged_DF(j,:));
                corr_temp(j)=temp(1,2);
            end
        end
        [max_corr idx_max]=max(corr_temp);
        if max_corr>0.75
            MergedIDX{counter}=[ToMerge{counter} idx_max];
            Merged_DF(i,:)=(Merged_DF(i,:)+Merged_DF(idx_max,:))/2;
            Merged_DF(idx_max,:)=[];
        end
    end
    counter=counter+1;
    progressbar(i/size(ALL_DF_02,1));
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
for i=1:8
    subplot(4,2,i);plot(Cmap(i,:))
    title(num2str(sum(idxKmeans==i)));
end

set(Fighandle, 'Position', [100, 100, 1200, 900]);
for i=1:8
    temp=find(idxKmeans==i);
    subplot(4,2,i);histogram(temp)
    title(num2str(sum(idxKmeans==i)));
end

coef=[Model_DF02(idx_015(idx)).coef];
coefficients_above015=[];
for idx=1:length(idx_015)
    coef=[Model_DF02(idx).coef];
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d)','tokens');temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
    for coef_idx=2:height(coef)
        if coef.pValue(coef_idx)<0.05
            coefficients_above015{idx,str2num(temp(coef_idx-1))}=coef.Estimate(coef_idx);
        end
    end
end
idxempty=cellfun('isempty',coefficients_above015);
coefficients_above015(idxempty)={0};
clearvars idxempty idx coef_idx coef
coefficients_above015=cell2mat(coefficients_above015);
[max_coef,idx_max]=max(coefficients_above015,[],2);

