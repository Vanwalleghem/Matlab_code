listOfROI = dir('*ROI_SpikePred.mat');
numberOfROI = numel(listOfROI);
x=load(listOfROI(1).name,'n_hat');
SpikeTraces=x.n_hat;

for j = 2:numberOfROI
	x=load(listOfROI(j).name,'n_hat');
    x=x.n_hat;
    SpikeTraces=cat(1,SpikeTraces,x);
end

listOfFluo = dir('*ROI_Raw.csv');
numberOfFluo = numel(listOfFluo);
x=csvread(listOfFluo(1).name);
FluoTraces=x';

for j = 2:numberOfFluo
	x=csvread(listOfFluo(j).name);
    FluoTraces=cat(1,FluoTraces,x');
end

x = linspace(1,size(H_ZS,2),size(H_ZS,2));
figure;
j=1;
for i = 1:3:size(H_ZS,1)
    if i==49
        subplot(size(H_ZS,1)/5,2,j);plot(x,H_ZS(i,:),x,H_ZS(i+1,:));
    else
        subplot(size(H_ZS,1)/5,2,j);plot(x,H_ZS(i,:),x,H_ZS(i+1,:),x,H_ZS(i+2,:));
    end
    j=j+1;
end

x = linspace(1,size(H_spike,2),size(H_spike,2));
figure;
j=1;
for i = 1:3:size(H_spike,1)
    if i==49
        subplot(size(H_spike,1)/5,2,j);plot(x,H_spike(i,:),x,H_spike(i+1,:));
    else
        subplot(size(H_spike,1)/5,2,j);plot(x,H_spike(i,:),x,H_spike(i+1,:),x,H_spike(i+2,:));
    end
    j=j+1;
end

x = linspace(1,size(H_ZS,2),size(H_ZS,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1049, 895]);
j=1;
for i = 1:3:size(mean_NMF_Fluo_sp,1)
    if i==49
        subplot(size(mean_NMF_Fluo_sp,1)/5,2,j);plot(x,mean_NMF_Fluo_sp(i,:),x,mean_NMF_Fluo_sp(i+1,:));
    else
        subplot(size(mean_NMF_Fluo_sp,1)/5,2,j);plot(x,mean_NMF_Fluo_sp(i,:),x,mean_NMF_Fluo_sp(i+1,:),x,mean_NMF_Fluo_sp(i+2,:));
    end    
    j=j+1;
end

figure;
subplot(3,1,1);
plot(x,mean_NMF_sp(6,:),x,mean_NMF_sp(8,:));
subplot(3,1,2);
plot(x,mean_NMF_sp(15,:),x,mean_NMF_sp(18,:));
subplot(3,1,3);
plot(x,mean_NMF_sp(19,:),x,mean_NMF_sp(21,:),x,mean_NMF_sp(36,:));

NMF_sp=[];
idx_sp=[];
for i = 1:size(H_spike,1)
    idx_sp{i}=find(W_spike(:,i)>1);
    NMF_sp{i}=Sp_infer(idx_sp{i},:);
end

mean_NMF_sp=zeros(size(H_spike,1),size(H_spike,2),'double');
for i = 1:size(H_spike,1)
    mean_NMF_sp(i,:)=mean(NMF_sp{i},1);
end

clear NMF_sp;

NMF_sp_Fluo=[];
idx_sp=[];
for i = 1:size(H_spike,1)
    idx_sp{i}=find(W_spike(:,i)>1);
    NMF_sp{i}=Sp_infer(idx_sp{i},:);
end

mean_NMF_sp=zeros(size(H_spike,1),size(H_spike,2),'double');
for i = 1:size(H_spike,1)
    mean_NMF_sp(i,:)=mean(NMF_sp{i},1);
end

NMF_ZS=[];
idx_ZS=[];
mean_NMF_ZS=zeros(size(H_ZS,1),size(H_ZS,2),'double');
for i = 1:size(H_spike,1)
    idx_ZS{i}=find(W_ZS(:,i)>4);
    NMF_ZS=FluoTraces(idx_ZS{i},:);
    mean_NMF_ZS(i,:)=mean(NMF_ZS,1);
end


mean_NMF_ZS=zeros(size(H_ZS,1),size(H_ZS,2),'double');
for i = 1:size(H_ZS,1)
    mean_NMF_ZS(i,:)=mean(NMF_ZS{i},1);
end

xcorr_ZS=zeros(size(mean_NMF_sp,1),size(mean_NMF_sp,1),'double');
for i = 1:size(mean_NMF_sp,1)
    for j = 1:size(mean_NMF_sp,1)
        if j>i
            t=corrcoef(mean_NMF_ZS(i,:)',mean_NMF_ZS(j,:)');
            xcorr_ZS(i,j)=t(2,1);
        else
            xcorr_ZS(i,j)=0;
        end    
    end    
end
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 1049, 895]);
imagesc(xcorr_ZS,[0 1]);

NMF_ZS=[];
idx_ZS=[];
mean_NMF_ZS=zeros(size(H_ZS,1),size(H_ZS,2),'double');
for i = 1:size(H_ZS,1)
    idx_ZS{i}=find(W_ZS(:,i)>2);
    NMF_ZS=ZS(idx_ZS{i},:);
    mean_NMF_ZS(i,:)=mean(NMF_ZS,1);   
end

mean_NMF_ZS=zeros(size(H_spike,1),size(H_spike,2),'double');
for i = 1:size(H_spike,1)
    mean_NMF_ZS(i,:)=mean(NMF_ZS{i},1);    
end

NMF_Fluo_sp=[];
for i = 1:size(H_spike,1)    
    NMF_Fluo_sp{i}=Fluo_select(idx_sp{i},:);
end
mean_NMF_Fluo_sp=zeros(size(H_spike,1),size(H_spike,2),'double');
for i = 1:size(H_spike,1)
    mean_NMF_Fluo_sp(i,:)=mean(NMF_Fluo_sp{i},1);    
end

xcorr_ZS=zeros(size(H_spike,1),size(H_spike,1),'double');
for i = 1:size(H_spike,1)
    for j = 1:size(H_spike,1)
        if j>i
            t=corrcoef(mean_NMF_ZS(i,:)',mean_NMF_ZS(j,:)');
            xcorr_ZS(i,j)=t(2,1);            
        else
            xcorr_ZS(i,j)=0;
        end    
    end    
end
figure;imagesc(xcorr_ZS,[0 1]);


figure;
j=1;
for i = 1:3:size(H_spike,1)
    subplot(size(H_spike,1)/5,2,j);plot(x,mean(NMF_sp{i}),x,mean(NMF_sp{i+1}),x,mean(NMF_sp{i+2}));
    j=j+1;
end


Results=cell(size(Fluo_select,1),6);
for i=1:size(Fluo_select,1)
    Results{i,1},Results{i,2},Results{i,3},Results{i,4},Results{i,5},Results{i,6}=constrained_foopsi(Fluo_select(i,:));
end


listOfMasks2 = dir('Maskb_*.tif');
numberOfMasks2 = numel(listOfMasks2);
Y2(numberOfMasks2,1)=0;
for i = 1:numberOfMasks2
	image=imread(listOfMasks2(i).name);
    A2{i,1}=image;
    Y2(i)=max(max(image));
end

parfor i=1:size(A,1)
    image=zeros(size(A{i,1}),'uint16');
    counter=1;
    idx=unique(A{i,1});
    for j=idx'        
        image(A{i,1}==j)=counter;
        counter=counter+1;
    end
    A2{i,1}=image;
end

Y2(numberOfMasks,1)=0;
parfor i=1:size(A2,1)
    Y2(i)=max(max(A2{i,1}));    
end


idx_linreg(size(temp,1),1)=0;
correction=0;
counter=1;
for i = 1:size(mylinreg,1)
    if size(mylinreg{i,3},1) > 0 
        for j = 1:size(mylinreg{i,3},1)
            idx_linreg(counter)=correction+mylinreg{i,3}(j);
            counter=counter+1;
        end
    end
    correction=correction+size(mylinreg{i,2},1);
end

Y3(numberOfMasks,1)=0;
parpool(4)
for i=1:size(Filelistb,2)
    name=strcat('Maskb_',Filelistb{1,i});
    image=imread(name);
    image2=zeros(size(image),'uint16');
    counter=1;
    idx=unique(image);
    for j=idx'        
        image2(image==j)=counter;
        counter=counter+1;
    end
    A3{i,1}=image2;
    Y3(i)=max(max(image2));
end

idx_select_compSpike=[];
for i = 1:numel(select_comps)
    idx_select_compSpike{i}=find(W_spike(:,select_comps(i))>1);
end

mean_select_compSpike=zeros(numel(select_comps),size(Fluo_select,2));
for i = 1:numel(select_comps)
    mean_select_compSpike(i,:)=mean(Fluo_select(idx_select_compSpike{i},:));
end

figure;
for i = 1:numel(select_comps)
    subplot(numel(select_comps),1,i);plot(x,mean_select_compSpike(i,:));
end

MI_selectComp=[];
for i = 1: numel(idx_select_compSpike)
    for j = 1: numel(idx_select_compSpike)
        test=intersect(idx_select_compSpike{i},idx_select_compSpike{j});
        MI_selectComp{i,j}=2*numel(test)/(numel(idx_select_compSpike{i})+numel(idx_select_compSpike{j}));
    end
end
MI_selectComp=cell2mat(MI_selectComp);

sum_explained(size(explained,1),1)=0;
sum_explained_sp(size(explained,1),1)=0;
sum_explained(1)=explained(1);
sum_explained_sp(1)=explained_sp(1);
for i = 2:size(explained,1)
    sum_explained(i)=sum_explained(i-1)+explained(i);
    sum_explained_sp(i)=sum_explained_sp(i-1)+explained_sp(i);
end