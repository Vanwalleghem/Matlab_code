load('GCaMP6f2-FluoTraces-Maskb.mat')
Nb_corr_CTRL=load('GCaMP6f2-FluoTraces-Maskb-SelectStrict.mat','Nb_corr');
Nb_corr_CTRL=Nb_corr_CTRL.Nb_corr;
idx_corr=find(Nb_corr_CTRL>5);
CTRL=FluoTraces(idx_corr);
CTRL=FluoTraces(idx_corr,:);CTRL=double(CTRL);
clearvars FluoTraces
load('neomycin-FluoTraces-Maskb.mat', 'FluoTraces')
load('neomycin-FluoTraces-Maskb-SelectCorr.mat', 'Nb_corr')
idx_Neo=find(Nb_corr>5);
Neo=double(FluoTraces(idx_Neo,:));
ALL_DF=[CTRL;Neo];
ALL_DF=DeltaF2(double(ALL_DF),21,11);
max_DF=max(ALL_DF,[],2);min_DF=min(ALL_DF,[],2);
idx_Threshold_5to200=find(max_DF>0.05 & max_DF<2 & min_DF>-0.1);
%Plane_correction=[0 0 10 0 10 0 0 0];%Careful that planes are converted to neg values so 20um becomes -20 so if 60um F5 = 50um F8 -> the correction = +10um
%Plane_correction_Neo=[-10 -20 -10 0 -30 -20 0];
%Plane_correction_Neo=[0 -10 0 10 20 0 20];
%Plane_correction=[-10 50 0 -10 20 10 -20 0];
Plane_correction=[0 40 10 0 10 10 -10 0];
Plane_correction_Neo=[0 -10 0 10 0 -10 10];

ALL_DF=ALL_DF(idx_Threshold_5to200,:);
clearvars Nb_corr FluoTraces

parfor i=1:size(ALL_DF,1)
mdl=stepwiselm(aud8freq',ALL_DF(i,:),'linear','Criterion','bic','Upper','linear','verbose',0);
Mode_ALLDF(i).coef=mdl.Coefficients;
Mode_ALLDF(i).MSE=mdl.MSE;
Mode_ALLDF(i).Fitted=mdl.Fitted;
Mode_ALLDF(i).rsquared=mdl.Rsquared.Adjusted;
end

rsq=[Mode_ALLDF.rsquared];
idx_015=find(rsq>0.15);
ALL_DF_02=ALL_DF(idx_015,:);



Merged_DF=ALL_DF_02;%   ALL_DF_02=Dataset (neurons x time)
MergedIDX=cell(1,length(Merged_DF)); % Where the IDs of merged ROIs will be stored


counter=1;
progressbar % can be skipped if you haven't progressbar installed    
while counter<=size(Merged_DF,1) % can also be for counter=1:size(Dataset,1)    
    if isnan(Merged_DF(counter,1)) % if the timeserie has been merged, it's set to nan, so this prevents spending time on it
        max_corr=0;
        MergedIDX{counter}=[];
    else
        max_corr=1;    % hack to get the equivalent of a do... while loop (at least one iteration)
        if isempty(MergedIDX{counter})
            MergedIDX{counter}=counter;
        end
    end    
    while max_corr>0.85 % standard threshold of Bianco et al
        corr_temp=zeros(size(Merged_DF,1),1);        
        parfor j=1:size(Merged_DF,1)    % parallelized, can replace by for loop if not needed
        	if j>counter && ~isnan(Merged_DF(j,1))
                temp=corrcoef(Merged_DF(counter,:), Merged_DF(j,:));
                corr_temp(j)=temp(1,2);                                
            end
        end        
        [max_corr idx_max]=nanmax(corr_temp);
        if max_corr>0.85
            if max_corr>0.95 % shortcut to merge in one go all time series with >0.85 correlation, can be skipped or changed
                idx_max=find(corr_temp>0.95);
                MergedIDX{counter}=[MergedIDX{counter} idx_max'];
            else
                MergedIDX{counter}=[MergedIDX{counter} idx_max]; % merge things one by one, so it's long
            end
            Merged_DF(counter,:)=nanmean(ALL_DF_02(MergedIDX{counter},:),1); %average from dataset, easier than adding one at a time with a factor                
            Merged_DF(idx_max,:)=nan;
            corr_temp(idx_max)=nan;            
        end
    end
    counter=counter+1;
    progressbar(counter/size(Merged_DF,1)); % can be skipped if you haven't progressbar installed    
end
clearvars i j idx_max counter max_corr temp

numMerged=zeros(numel(MergedIDX),1);
for i=1:numel(MergedIDX)
    numMerged(i)=numel(MergedIDX{i});
end
%figure;histogram(numMerged(numMerged>10));

temp=find(numMerged>200);
MergedComps=[];
counter=1;
for i=temp'
    MergedComps(counter).mean=mean(ALL_DF_02(MergedIDX{i},:)*100,1);
    MergedComps(counter).STD=std(ALL_DF_02(MergedIDX{i},:)*100,1,1);
    MergedComps(counter).NB=ones(1,size(ALL_DF_02,2))*numel(MergedIDX{i});
    MergedComps(counter).Distrib=idx_fish_ALL_02(MergedIDX{i});
    counter=counter+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1000]);
counter=1;xplot=floor(sqrt(length(temp)));yplot=ceil(length(temp)/xplot);
for i=1:length(temp)
    subplot(xplot,yplot,i);plot(MergedComps(i).mean);
    title(num2str(MergedComps(i).NB(1)));    
end
clearvars temp i


Test=zeros(length(temp),700);
counter=1;
for i=temp'
    Test(counter,:)=mean(ALL_DF_02(MergedIDXc{i},:)*100,1);
    mdl=stepwiselm(aud8freq',Test(counter,:),'linear','verbose',0);
    model_clusters(counter).coef=mdl.Coefficients;
    model_clusters(counter).MSE=mdl.MSE;
    model_clusters(counter).Fitted=mdl.Fitted;
    model_clusters(counter).rsquared=mdl.Rsquared.Adjusted;
    counter=counter+1;
end
clearvars counter i counter

PrismTrans=zeros(3*numel(temp),size(ALL_DF_02,2)); % to transfer data in prism as mean, std and n
counter=1;
for i=1:numel(temp)
    PrismTrans(counter,:)=MergedComps(i).mean;
    PrismTrans(counter+1,:)=MergedComps(i).STD;
    PrismTrans(counter+2,:)=MergedComps(i).NB;
    counter=counter+3;
end

counter=1;
idx_fish_CTRL=zeros(length(idx_corr),1);
[Fish,~]=regexp(files_info_CTRL.names(1),'F(\d+)','tokens','match');Fish=str2num(Fish{1}{1}{1});
idx_fish_CTRL(1:files_info_CTRL.index(1))=Fish;
for i=2:length(files_info_CTRL.index)
	[Fish,~]=regexp(files_info_CTRL.names(i),'F(\d+)','tokens','match');Fish=str2num(Fish{1}{1}{1});	
    idx_fish_CTRL(files_info_CTRL.index(i-1):files_info_CTRL.index(i))=Fish;
end

idx_fish_Neo=zeros(length(idx_corr_bool),1);
[Fish,~]=regexp(files_info_Neo.names(1),'F(\d+)','end');Fish=Fish{1,1};
tempcomp=files_info_Neo.names(1);tempcomp=tempcomp{1,1}(1:Fish);
counter2=1;
idx_fish_Neo(1:files_info_Neo.index(1))=counter2;
for i=2:length(files_info_Neo.index)
	[Fish,~]=regexp(files_info_Neo.names(i),'F(\d+)','end');Fish=Fish{1,1};
	tempcomp2=files_info_Neo.names(i);tempcomp2=tempcomp2{1,1}(1:Fish);      
	if ~strcmp(tempcomp,tempcomp2)
        tempcomp=tempcomp2;
        counter2=counter2+1;
	end
    idx_fish_Neo(files_info_Neo.index(i-1):files_info_Neo.index(i))=counter2;
end
figure;histogram(idx_fish_Neo);

idx_fish_CTRL_corr=idx_fish_CTRL(idx_corr);
idx_fish_Neo_corr=idx_fish_Neo(idx_corr_bool);
idx_fish_Neo_corr=idx_fish_Neo_corr(idx_Threshold_5to200_bool);
idx_fish_ALL=[idx_fish_CTRL_corr ; idx_fish_Neo_corr+max(idx_fish_CTRL_corr)];
idx_fish_ALL_02=idx_fish_ALL(idx_02);
idx_fish_ALL_015=idx_fish_ALL(idx_015);


Centroids_CTRL_select=Centroids_CTRL(idx_corr_CTRL,:);
Centroids_Neo_select=Centroids_Neo(idx_corr_Neo,:);
Centroids_ALL=[Centroids_CTRL_select;Centroids_Neo_select];
Centroids_ALL=Centroids_ALL(idx_Threshold_5to200,:);
Centroids_015=Centroids_ALL(idx_015,:);

numMerged=zeros(numel(MergedIDX),1);
for i=1:numel(MergedIDX)
    numMerged(i)=numel(MergedIDX{i});
end
temp=find(numMerged>50);
Representative_clust=[];
for i=1:numel(temp)
    temp_idx=MergedIDX{temp(i)};
    temp_fish=idx_fish_ALL_02(temp_idx);
    if numel(unique(temp_fish))>=5   
        nbcells=[];
        for j=1:max(idx_fish_ALL_02)
            nbcells{j}=sum(temp_fish==j);
        end
        nbcells=cell2mat(nbcells);
        nbcells=sum(nbcells>=5);
        if nbcells>=5
            Representative_clust=[Representative_clust temp(i)];
        end
    end
end
clearvars temp i temp_idx temp_fish nbcells

PrismTrans=zeros(3*numel(Representative_clust),size(ALL_DF_02,2)); % to transfer data in prism as mean, std and n
counter=1;
for i=1:numel(Representative_clust)
    PrismTrans(counter,:)=mean(ALL_DF_02(MergedIDX{Representative_clust(i)},:)*100,1);
    PrismTrans(counter+1,:)=std(ALL_DF_02(MergedIDX{Representative_clust(i)},:)*100,1,1);
    PrismTrans(counter+2,:)=ones(1,size(ALL_DF_02,2))*numel(MergedIDX{Representative_clust(i)});
    counter=counter+3;
end

PrismTransDistrib=zeros(numel(Representative_clust),15);
for i=1:numel(Representative_clust)
    temp=idx_fish_ALL_02(MergedIDX{Representative_clust(i)});
    temp2=zeros(1,15);
    for j=1:15
        temp2(j)=sum(temp==j);
    end
    temp2=temp2./sum(temp2);
    PrismTransDistrib(i,:)=temp2;
end
PrismTransDistrib=PrismTransDistrib*100;
clearvars i j temp temp2

MeanClust=zeros(length(Representative_clust),700);
counter=1;
for idx=Representative_clust
MeanClust(counter,:)=mean(ALL_DF_02(MergedIDX{idx},:),1);
counter=counter+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1000]);
counter=1;xplot=floor(sqrt(length(Representative_clust)));yplot=ceil(length(Representative_clust)/xplot);
for i=1:length(Representative_clust)
    subplot(xplot,yplot,i);plot(MeanClust(i,:));
    title(num2str(numel(ALL_DF_02(MergedIDX{Representative_clust(i)}))));    
end


%For second round
Merged_DFb=Merged_DF;
MergedIDXb=MergedIDX;

progressbar
counter=length(Representative_clust);
while counter>=1
    if isnan(Merged_DF(Representative_clust(counter),1)) % if the timeserie has been merged, it's set to nan, so this prevents spending time on it
        max_corr=0;
        MergedIDX{Representative_clust(counter)}=[];
    else
        max_corr=1;    % hack to get the equivalent of a do... while loop (at least one iteration)
        if isempty(MergedIDX{Representative_clust(counter)})
            MergedIDX{Representative_clust(counter)}=Representative_clust(counter);
        end
    end    
    while max_corr>0.85 % standard threshold of Bianco et al
        corr_temp=zeros(size(Merged_DF,1),1);        
        parfor j=1:size(Merged_DF,1)    % parallelized, can replace by for loop if not needed
        	if ~any(ismember(j,Representative_clust)) && ~isnan(Merged_DF(j,1))
                temp=corrcoef(Merged_DF(Representative_clust(counter),:), Merged_DF(j,:));
                corr_temp(j)=temp(1,2);                                
            end
        end        
        [max_corr idx_max]=nanmax(corr_temp);
        if max_corr>0.85
            if max_corr>0.95 % shortcut to merge in one go all time series with >0.85 correlation, can be skipped or changed
                idx_max=find(corr_temp>0.95);
                for idx=1:length(idx_max)
                    if isempty(MergedIDX{idx_max(idx)})
                        MergedIDX{Representative_clust(counter)}=[MergedIDX{Representative_clust(counter)} idx_max(idx)];
                        MergedIDX{idx_max(idx)}=[];
                    else
                        MergedIDX{Representative_clust(counter)}=[MergedIDX{Representative_clust(counter)} MergedIDX{idx_max(idx)}];
                        MergedIDX{idx_max(idx)}=[];
                    end
                end
            else
                if isempty(MergedIDX{idx_max})
                    MergedIDX{Representative_clust(counter)}=[MergedIDX{Representative_clust(counter)} idx_max]; % merge things one by one, so it's long
                    MergedIDX{idx_max}=[];
                else
                    MergedIDX{Representative_clust(counter)}=[MergedIDX{Representative_clust(counter)} MergedIDX{idx_max}];
                    MergedIDX{idx_max}=[];
                end
            end
            Merged_DF(Representative_clust(counter),:)=nanmean(ALL_DF_02(MergedIDX{Representative_clust(counter)},:),1); %average from dataset, easier than adding one at a time with a factor                
            Merged_DF(idx_max,:)=nan;
            corr_temp(idx_max)=nan;            
        end
    end
    counter=counter-1;
    progressbar((length(Representative_clust)-counter)/length(Representative_clust)); % can be skipped if you haven't progressbar installed    
end
clearvars i j idx_max counter max_corr temp

idx_GoodBetasSelect_CTRL=[];temp2=[];
for i=1:length(Representative_clust)
    temp=MergedIDX{Representative_clust(i)};
    idx_GoodBetasSelect_CTRL=[idx_GoodBetasSelect_CTRL temp];
    temp2=[temp2 ones(size(temp))*Representative_clust(i)];
end
idx_GoodBetasSelect_CTRL=[idx_GoodBetasSelect_CTRL; temp2];
idx_GoodBetasSelect_CTRL=idx_GoodBetasSelect_CTRL';

idx_fish_RepClusters=[];
for i=1:numel(Representative_clust)
    idx_fish_RepClusters=[idx_fish_RepClusters; idx_fish_ALL_02(MergedIDX{Representative_clust(i)})];
end

idx_rep_clust=[];
for i=1:numel(Representative_clust)
    idx_rep_clust=[idx_rep_clust MergedIDX{Representative_clust(i)}];
end

% angles=degtorad([-3.4 -14.3 -1.5 -1.6 92 98 -2 -2.5]);
% Plane_correction=[-10 50 0 -10 20 10 -20 0];
% %Plane_correction=[0 0 0 0 0 0 0 0];
% Centroids_CTRL_corrected={};
% for Fish=1:8
%     idx=find(idx_fish_RepClusters==Fish);%idx=idx_01(idx);    
%     temp_x=double(Centroids_015(idx,2));
%     temp_y=double(Centroids_015(idx,1));
%     Centroids_CTRL_corrected{Fish,5}=idx_GoodBetasSelect_CTRL(idx,2);    
%     temp_planes=idx_Plane_02(idx);list_planes=unique(temp_planes);
%     for Plane=list_planes'
%         idx_plane=find(temp_planes==Plane);
%         index = find(ismember(TransfoMatricescleared(:,1), strcat(num2str(Fish),'_',num2str(Plane))));
%         if index
%             index_plane=(Plane-40)/10;
%             if index_plane<1
%                 index_plane=1;
%             end
%             temp_x(idx_plane) = temp_x(idx_plane)+TransfoMatricescleared{index,2};
%             temp_y(idx_plane) = temp_y(idx_plane)+TransfoMatricescleared{index,3};
%         end
%     end
%     temp_x=(temp_x-prctile(temp_x,2))/prctile(temp_x,98);
%     temp_y=(temp_y-prctile(temp_y,2))/prctile(temp_y,98);
% 	temp_xb=temp_x*cos(angles(Fish))-temp_y*sin(angles(Fish));temp_xb=temp_xb+abs(min(temp_xb));
% 	temp_yb=temp_y*cos(angles(Fish))+temp_x*sin(angles(Fish));temp_yb=temp_yb+abs(min(temp_yb));    
%     Centroids_CTRL_corrected{Fish,1}=temp_xb;
%     Centroids_CTRL_corrected{Fish,2}=temp_yb;        
%     Centroids_CTRL_corrected{Fish,3}=temp_planes-Plane_correction(Fish);
%     Centroids_CTRL_corrected{Fish,4}=ones(size(temp_xb))*Fish;    
% end
% clearvars idx counter temp_x temp_y temp_xb temp_yb Fish Plane i index index_plane j Plane_correction angles idx_plane temp_Fish temp_planes temp
% X_coord=vertcat(Centroids_CTRL_corrected{:,1});
% Y_coord=vertcat(Centroids_CTRL_corrected{:,2});
% Z_coord=vertcat(Centroids_CTRL_corrected{:,3});
% XYZ_coord_CTRL=horzcat(X_coord,Y_coord,Z_coord); clearvars X_coord Y_coord Z_coord
% Fish_index=vertcat(Centroids_CTRL_corrected{:,4});
% coef_index=vertcat(Centroids_CTRL_corrected{:,5});

angles=degtorad([-3.4 -14.3 -1.5 -1.6 92 98 -2 -2.5]);
%Plane_correction=[-10 50 0 -10 20 10 -20 0];
Centroids_CTRL_corrected=zeros(size(Centroids_CTRL,1),4);
for i=1:length(files_info_CTRL.index)
%for i=[1:1:93 141:1:190]
%for i=[1:1:72 141:1:190]
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
    index = find(ismember(TransfoMatricescleared(:,1), strcat(num2str(Fish),'_',num2str(Plane))));%index_plane=(Plane-40)/10;
%     if index_plane<1
%     	index_plane=1;
%     end    
    temp_x = temp_x+TransfoMatricescleared{index,2};%+CenteringMatrix(Fish,index_plane,1);
    temp_y = temp_y+TransfoMatricescleared{index,3};%+CenteringMatrix(Fish,index_plane,2);
	temp_xb=temp_x*cos(angles(Fish))-temp_y*sin(angles(Fish));%temp_xb=temp_xb+abs(min(temp_xb));
	temp_yb=temp_y*cos(angles(Fish))+temp_x*sin(angles(Fish));%temp_yb=temp_yb+abs(min(temp_yb));
    %temp_xb=double(temp_xb);temp_yb=double(temp_yb);
%     if min(temp_xb)<0
%         temp_xb=temp_xb+1080;
%     end
%     if min(temp_yb)<0
%         temp_yb=temp_yb+1280;
%     end
    %temp_xb=(temp_xb-prctile(temp_xb,2))/prctile(temp_xb,98);
    %temp_yb=(temp_yb-prctile(temp_yb,2))/prctile(temp_yb,98);
    Centroids_CTRL_corrected(counter:idx,1)=temp_xb;
    Centroids_CTRL_corrected(counter:idx,2)=temp_yb;        
    Centroids_CTRL_corrected(counter:idx,3)=-Plane;
    Centroids_CTRL_corrected(counter:idx,4)=Fish;
    counter=idx+1;
end
clearvars idx counter temp_x temp_y Fish Plane i

Centroids_CTRL_correctedb=zeros(size(Centroids_CTRL,1),4);
for Fish=1:8;
    idx=find(Centroids_CTRL_corrected(:,4)==Fish);
    temp_x=double(Centroids_CTRL_corrected(idx,1));temp_y=double(Centroids_CTRL_corrected(idx,2));temp_x=temp_x-min(temp_x);temp_y=temp_y-min(temp_y);
    temp_x=(temp_x-prctile(temp_x,2))/prctile(temp_x,98);
    temp_y=(temp_y-prctile(temp_y,2))/prctile(temp_y,98);
    if Fish<8
        [temp_x, temp_y]=transformPointsForward(Geometric_trans_CTRL_scaled{Fish},temp_y,temp_x);
    else
        temp=temp_x;temp_x=temp_y;temp_y=temp;
    end
    Centroids_CTRL_correctedb(idx,:)=[temp_x temp_y Centroids_CTRL_corrected(idx,3)+Plane_correction(Fish) Centroids_CTRL_corrected(idx,4)];
end
figure;scatter(Centroids_CTRL_correctedb(:,1),Centroids_CTRL_correctedb(:,2));
clearvars Fish temp_x temp_y temp

XYZ_coord_CTRL=Centroids_CTRL_correctedb(:,1:3);
Fish_index_CTRL=Centroids_CTRL_correctedb(:,4);
% coef_index=vertcat(Centroids_CTRL_corrected{:,5});

angles_neo=degtorad([10 7 0.69 -7.17 -2 6 -5.5]);
%Plane_correction_Neo=[0 -10 0 10 20 0 20];
filename=files_info_Neo.names(1);filename=filename{1};
[Fish,~]=regexp(filename,'F(\d+)','end');
tempcomp=filename(1:Fish);counter2=1;
Centroids_Neo_corrected=cell(7,4);%counter3=1;
for i=1:length(files_info_Neo.names)
%figure;counter2=4;
%for i=[60:1:80]
    %filename=neoTransfoMatrix{i,1};
    filename=files_info_Neo.names(i);filename=filename{1};
    test=find(cellfun( @(x) strcmp(x,filename), {neoTransfoMatrix{:,1}} ));
    if i==1
        counter=1;
    else
        counter=files_info_Neo.index(i-1)+1;
    end
	[Fish,~]=regexp(filename,'F(\d+)','end');
    tempcomp2=filename;tempcomp2=tempcomp2(1:Fish);
    if ~strcmp(tempcomp,tempcomp2)
        tempcomp=tempcomp2;
        counter2=counter2+1;
    end
    [Plane,~]=regexp(filename,'(\d+)um','tokens','match');Plane=str2num(Plane{1}{1});
    idx=files_info_Neo.index(i);  
    temp_x=Centroids_Neo(counter:idx,2);
    temp_y=Centroids_Neo(counter:idx,1);
    temp_x = temp_x+neoTransfoMatrix{test,2};%+CenteringMatrix(Fish,index_plane,1);
    temp_y = temp_y+neoTransfoMatrix{test,3};%+CenteringMatrix(Fish,index_plane,2);
%     temp_x=(temp_x-prctile(temp_x,2))/prctile(temp_x,98);
%     temp_y=(temp_y-prctile(temp_y,2))/prctile(temp_y,98);
	temp_xb=temp_x*cos(angles_neo(counter2))-temp_y*sin(angles_neo(counter2));%temp_xb=temp_xb+abs(min(temp_xb));
	temp_yb=temp_y*cos(angles_neo(counter2))+temp_x*sin(angles_neo(counter2));%temp_yb=temp_yb+abs(min(temp_yb));
%      temp_xb=(temp_xb-prctile(temp_xb,2))/prctile(temp_xb,98);
%      temp_yb=(temp_yb-prctile(temp_yb,2))/prctile(temp_yb,98);
    temp_planes=-idx_Plane_Neo(counter:idx);
    Centroids_Neo_corrected{counter2,1}=[Centroids_Neo_corrected{counter2,1}; temp_xb];
    Centroids_Neo_corrected{counter2,2}=[Centroids_Neo_corrected{counter2,2}; temp_yb];        
    Centroids_Neo_corrected{counter2,3}=[Centroids_Neo_corrected{counter2,3}; temp_planes];
    Centroids_Neo_corrected{counter2,4}=[Centroids_Neo_corrected{counter2,4}; ones(size(temp_xb))*counter2];    
    %counter=idx+1;
    %subplot(5,5,counter3);scatter(temp_xb,temp_yb);
    %counter3=counter3+1;
end
clearvars idx counter temp_x temp_y Fish Plane i idx_plane i counter2 temp_xb temp_yb temp_planes test tempcomp tempcomp2 C ia ib test

Centroids_Neo_correctedb=zeros(size(Centroids_Neo,1),4);
for Fish=1:7;
    idx=find(Fish_index_Neo==Fish);
    temp_x=double(Centroids_Neo_corrected{Fish,1});temp_y=double(Centroids_Neo_corrected{Fish,2});temp_x=temp_x-min(temp_x);temp_y=temp_y-min(temp_y);
    temp_x=(temp_x-prctile(temp_x,2))/prctile(temp_x,98);
    temp_y=(temp_y-prctile(temp_y,2))/prctile(temp_y,98);
    [temp_x, temp_y]=transformPointsForward(Geometric_trans_Neo_scaled{Fish},temp_y,temp_x);   
    Centroids_Neo_correctedb(idx,:)=[temp_x temp_y Centroids_Neo_corrected{Fish,3}+Plane_correction_Neo(Fish) Centroids_Neo_corrected{Fish,4}];
end
%figure;scatter(Centroids_Neo_correctedb(:,1),Centroids_Neo_correctedb(:,2));
%figure;scatter(Centroids_CTRL_correctedb(:,1),Centroids_CTRL_correctedb(:,2));

% X_coord=vertcat(Centroids_Neo_corrected{:,1});
% Y_coord=vertcat(Centroids_Neo_corrected{:,2});
% Z_coord=vertcat(Centroids_Neo_corrected{:,3});
% XYZ_coord_Neo=horzcat(X_coord,Y_coord,Z_coord); clearvars X_coord Y_coord Z_coord
Fish_index_Neo=Centroids_Neo_correctedb(:,4);
XYZ_coord_Neo=Centroids_Neo_correctedb(:,1:3);
% 
% Centroids_CTRL_select=Centroids_CTRL(idx_corr_CTRL,:);Centroids_CTRL_select=[Centroids_CTRL_select zeros(length(Centroids_CTRL_select),1)];
% Centroids_Neo_select=XYZ_coord_Neo(idx_corr_Neo,:);
% Centroids_ALL=[double(Centroids_CTRL_select);Centroids_Neo_select];
% Centroids_015=Centroids_ALL(idx_015,:);

Centroids_CTRL_select=XYZ_coord_CTRL(idx_corr,:);
Centroids_Neo_select=XYZ_coord_Neo(idx_Neo,:);
ALL_Cent_select=[Centroids_CTRL_select;Centroids_Neo_select];
ALL_Cent_select=ALL_Cent_select(idx_Threshold_5to200,:);
Centroids_015=ALL_Cent_select(idx_015,:);

Fish_index_CTRL_select=Fish_index_CTRL(idx_corr);
Fish_index_Neo_select=Fish_index_Neo(idx_Neo);
Fish_index_ALL=[Fish_index_CTRL_select; Fish_index_Neo_select+max(Fish_index_CTRL_select)];
Fish_index_ALL=Fish_index_ALL(idx_Threshold_5to200);
Fish_index_015=Fish_index_ALL(idx_015);

% XYZ_coord_Neo_select={};
% for i=1:length(Representative_clust)
%     temp=Centroids_015(MergedIDX{Representative_clust(i)},:);
%     temp2=Fish_index_015(MergedIDX{Representative_clust(i)});
%     temp(temp(:,3)==0,:)=[];
%     temp2(temp2==0)=[];
%     XYZ_coord_Neo_select{i,1}=temp;
%     XYZ_coord_Neo_select{i,2}=temp2;
% end
% 
XYZ_coord_select={};
Fish_index_select={};
for i=1:length(Representative_clust)
    XYZ_coord_select{i}=Centroids_015(MergedIDX{Representative_clust(i)},:);
    Fish_index_select{i}=Fish_index_015(MergedIDX{Representative_clust(i)});
end

Ganglia_Limits_CTRL=[0.5 0.62;0.43 0.64;0.45 0.45;0.47 0.6;0.45 0.57;0.53 0.7;0.5 0.7;0.48 0.6];
Ganglia_Limits_Neo=[0.6 0.5; 0.5 0.5; 0.6 0.6; 0.4 0.3;0.5  0.3;0.4 0.4; 0.55 0.47];
Ganglia_Limits=[Ganglia_Limits_CTRL; Ganglia_Limits_Neo];
ganglia_coord_select=cell(1,length(Representative_clust));
Fish_index_gang_select=cell(1,length(Representative_clust));
for i=1:length(Representative_clust)
    for Fish=1:15
        idx_fish=find(Fish_index_select{i}==Fish);
        ganglia_coord_select{i}=[ganglia_coord_select{i}; XYZ_coord_select{i}(XYZ_coord_select{i}(idx_fish,1)>Ganglia_Limits(Fish,1) & XYZ_coord_select{i}(idx_fish,2)>Ganglia_Limits(Fish,2) & XYZ_coord_select{i}(idx_fish,2)<Ganglia_Limits(Fish,2)+0.4 & XYZ_coord_select{i}(idx_fish,3)<200,:)];
        Fish_index_gang_select{i}=[Fish_index_gang_select{i}; Fish_index_select{i}(XYZ_coord_select{i}(idx_fish,1)>Ganglia_Limits(Fish,1) & XYZ_coord_select{i}(idx_fish,2)>Ganglia_Limits(Fish,2) & XYZ_coord_select{i}(idx_fish,2)<Ganglia_Limits(Fish,2)+0.4 & XYZ_coord_select{i}(idx_fish,3)<200)];
    end
end
clearvars Ganglia_Limits_CTRL Ganglia_Limits_Neo i Fish

angles_ganglia_CTRL=degtorad([-52 -52 -55 -49 -45 -46 -41 -40]);
angles_ganglia_Neo=degtorad([-44 -50 -48 -37 -64 -49 -41]);
angles_ganglia=[angles_ganglia_CTRL angles_ganglia_Neo];
ganglia_coord_corrected_select={};
for Fish=1:15
    for j=1:length(Representative_clust)
        idx=find(Fish_index_gang_select{j}==Fish);
        temp_y=ganglia_coord_select{j}(idx,2)-Ganglia_Limits(Fish,2);temp_z=double(ganglia_coord_select{j}(idx,3))/200;
        ganglia_coord_corrected_select{Fish,j,2}=cos(angles_ganglia(Fish))*temp_y+sin(angles_ganglia(Fish))*temp_z;
        ganglia_coord_corrected_select{Fish,j,3}=cos(angles_ganglia(Fish))*temp_z-sin(angles_ganglia(Fish))*temp_y;
        ganglia_coord_corrected_select{Fish,j,1}=ganglia_coord_select{j}(idx,1);        
    end
end
clearvars temp Fish temp_y temp_z angles_ganglia_CTRL angles_ganglia_Neo j idx Ganglia_Limits filename

i=cellfun(@length,XYZ_coord_select);i=sum(i);
min_x=Centroids_Neo_corrected{:,1};max_x=max(min_x);min_x=min(min_x);
min_y=Centroids_Neo_corrected{:,2};max_y=max(min_y);min_y=min(min_y);
temp=zeros(i,4);counter=0;
for i=1:length(Representative_clust)
    temp(1+counter:counter+length(XYZ_coord_select{i}),1:3)=XYZ_coord_select{i};
    temp(1+counter:counter+length(XYZ_coord_select{i}),4)=ones(length(XYZ_coord_select{i}),1)*i;
    counter=counter+length(XYZ_coord_select{i});
end
temp(:,1)=0.6452*(max_x-min_x)*temp(:,1);
temp(:,2)=0.6452*(max_y-min_y)*temp(:,2);
clearvars i min_x min_y counter temp2 temp3 temp_xb temp_yb



PrismGangliaTrans=zeros(1,3*length(Representative_clust));counter=0;
for i=1:length(Representative_clust)
    PrismGangliaTrans(counter+1)=mean(vertcat(ganglia_coord_corrected_select{:,i,3}));
    PrismGangliaTrans(counter+2)=std(vertcat(ganglia_coord_corrected_select{:,i,3}));
    %PrismGangliaTrans(counter+3)=numel(vertcat(ganglia_coord_corrected_select{:,i,2}));    
    PrismGangliaTrans(counter+3)=numel(unique(Fish_index_gang_select{i}));
    counter=counter+3;
end
clearvars i counter

PrismTranscoord=zeros(1,3*length(Representative_clust));counter=0;
for i=1:length(Representative_clust)
    PrismTranscoord(counter+1)=mean(XYZ_coord_select{i}(:,3));
    PrismTranscoord(counter+2)=std(XYZ_coord_select{i}(:,3));    
    PrismTranscoord(counter+3)=numel(unique(Fish_index_select{i}));
    counter=counter+3;
end
clearvars i counter

Z_coord_corrected_select={};
for Fish=1:15
    for j=1:length(Representative_clust)
        idx=find(Fish_index_select{j}==Fish);
        temp_y=XYZ_coord_select{j}(idx,2);temp_z=double(XYZ_coord_select{j}(idx,3))/200;        
        Z_coord_corrected_select{Fish,j}=cos(angles_ganglia(Fish))*temp_z-sin(angles_ganglia(Fish))*temp_y;        
    end
end
clearvars i j Fish counter counter2 counter3 corr_temp

PrismTranscoef=zeros(1,3*8);counter=0;
for i=1:8
    idx=find(idx_max==i);
    PrismTranscoef(counter+1)=mean(Centroids_015(idx,3));
    PrismTranscoef(counter+2)=std(Centroids_015(idx,3));    
    PrismTranscoef(counter+3)=numel(unique(Fish_index_015(idx)));
    counter=counter+3;
end
clearvars i counter



Torus_Limits_CTRL=[0.5 0.62;0.43 0.64;0.45 0.45;0.47 0.6;0.45 0.57;0.53 0.7;0.5 0.7;0.48 0.6];
Torus_Limits_Neo=[0.6 0.5; 0.5 0.5; 0.6 0.6; 0.4 0.3;0.5  0.3;0.4 0.4; 0.55 0.47];
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

idx_final_CTRLb=[1:1:length(idx_corr_CTRL)];
idx_final_CTRLb=idx_final_CTRLb(idx_corr_CTRL);
idx_final_Neob=[1:1:length(idx_corr_bool)];
idx_final_Neob=idx_final_Neob(idx_corr_Neo);
idx_final_ALL=[idx_final_CTRLb idx_final_Neob; ones(1,length(idx_final_CTRLb)) ones(1,length(idx_final_Neob))*2];
idx_final_ALL=idx_final_ALL(:,idx_Threshold_5to200);
idx_final_015=idx_final_ALL(:,idx_015);
clearvars idx_final_CTRLb idx_final_Neob
idx_final_CTRL=zeros(length(idx_Plane_CTRL),1);
idx_final_Neo=zeros(length(idx_corr_bool),1);

for i=1:length(Representative_clust)
    idx_temp=MergedIDX{Representative_clust(i)};
    idx_temp=idx_final_015(:,idx_temp);
    idx_final_CTRL(idx_temp(1,find(idx_temp(2,:)==1)))=i;
    idx_final_Neo(idx_temp(1,find(idx_temp(2,:)==2)))=i;
end

figure;
for i=1:length(Representative_clust)
    subplot(length(Representative_clust),1,i);plot(MeanClust(i,:),'Color',colors(i,:));
end

colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
colors = colors*256;
for i=1:length(CTRL_Mask)
    B=CTRL_Mask{i};B=double(B);B(B==0)=NaN;
    filename=files_info_CTRL.names(i);
    filename=strcat('D:\Pictures\processed\Tonotropy\AVG_Cropped_CTRL\AVG_',filename{1});
    image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(max((image4)));image4=image4*128;image4=double(repmat(image4,1,1,3));image4=uint8(image4);
    %image=zeros(size(B,1),size(B,2),3);
    for j=min(min(B)):max(max(B))
        if any(idx_final_CTRL(j)>0)
            clust=idx_final_CTRL(j);
            for k=1:3              
              image2=image4(:,:,k);
              image2(B==j)=colors(clust,k);
              image4(:,:,k)=image2;
            end
        end
    end
    %image=uint8(image);
    filename=files_info_CTRL.names(i);
    name=strcat('Bianco_',filename{1});
    %filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\AVG_',sorted_files{i});
    %image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(prctile(image4,80));image4=image4*200;
    %image4=double(repmat(image4,1,1,3));image4=uint8(image4)+image;
    image4=uint8(image4);
    imwrite(image4,name,'tif');
end
clearvars image image2 image4 i j k name filename B clust;

for i=1:length(Neo_Mask)
    B=Neo_Mask{i};B=double(B);B(B==0)=NaN;
    filename=files_info_Neo.names(i);
    filename=strcat('D:\Pictures\processed\Tonotropy\neomycin\AVG_',filename{1});
    image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(max((image4)));image4=image4*128;image4=double(repmat(image4,1,1,3));image4=uint8(image4);
    %image=zeros(size(B,1),size(B,2),3);
    for j=min(min(B)):max(max(B))
        if any(idx_final_Neo(j)>0)
            clust=idx_final_Neo(j);
            for k=1:3              
              image2=image4(:,:,k);
              image2(B==j)=colors(clust,k);
              image4(:,:,k)=image2;
            end
        end
    end
    %image=uint8(image);
    filename=files_info_Neo.names(i);
    name=strcat('Bianco_',filename{1});
    %filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\AVG_',sorted_files{i});
    %image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(prctile(image4,80));image4=image4*200;
    %image4=double(repmat(image4,1,1,3));image4=uint8(image4)+image;
    image4=uint8(image4);
    imwrite(image4,name,'tif');
end
clearvars image image2 image4 i j k name filename B clust;
clearvars ans angle i j k counter counter2 Fish Plane

MeanClustb=MeanClust;
MeanClust=MeanClustb;
SelectRep=[3 2 5 1 4];

RRegionsTrans=zeros(length(Representative_clust)*4,700);counter=0;
for i=1:length(Representative_clust)
    temp=[gang_DF{1,SelectRep(i)};gang_DF{2,SelectRep(i)}];
    RRegionsTrans(counter+1,:)=mean(temp,1);
    temp=Hind_DF{SelectRep(i)};
    RRegionsTrans(counter+2,:)=mean(temp,1);
    temp=[Torus_DF{1,SelectRep(i)};Torus_DF{2,SelectRep(i)}];
    RRegionsTrans(counter+3,:)=mean(temp,1);
    temp=Thalamus_DF{SelectRep(i)};
    RRegionsTrans(counter+4,:)=mean(temp,1);
%     temp=Rest_DF{i};
%     RRegionsTrans(counter+5,:)=mean(temp,1);
%    counter=counter+5;    
    counter=counter+4;    
end

MeanClust=RRegionsTrans2;

lags=15;
%lag_cluster=zeros(((size(MeanClust,1)-1)*(size(MeanClust,1)))/2,1);
lag_cluster=zeros(size(MeanClust,1),size(MeanClust,1));
counter=1;
for i=1:size(MeanClust,1)
    for j=1:size(MeanClust,1)
         %if j>i
            [r,lag] = xcorr(MeanClust(j,:),MeanClust(i,:),lags,'coeff');
            f=fit(lag',r','gauss1');
            %lag_cluster(counter)=f.b1;
            lag_cluster(i,j)=f.b1;
            counter=counter+1;
        %end
    end
end
figure;bar(sum(lag_cluster,2))
%lag_tri=triu(lag_cluster);
figure;imagesc(lag_cluster,[-1 1]);colormap jet
