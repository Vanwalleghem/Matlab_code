counter=1;
for i = 1:size(files,1)    
    name=strcat('Maskb_',files{i,1});
    ispresent = cellfun(@(s) ~isempty(strfind(name,s)),listOfMasks);
    ispresent=find(ispresent(1,:)==1);
    image=zeros(size(A2{ispresent,1},1),size(A2{ispresent,1},2),3);
    image2=zeros(size(A2{ispresent,1,1}));
    image3=zeros(size(A2{ispresent,1,1}));
    filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\',name(7:end));
    image4=imread(filename);
    image4=single(image4);image4(image4<500)=0;image4=image4/max(max(image4));image4=image4*256;
    for j = 1:max(max(A2{ispresent,1}))
        if counter < index(i,1)
            image2(A2{ispresent,1,1}==j)=W_ZS_29(counter,1)*256;
            image3(A2{ispresent,1,1}==j)=W_ZS_30(counter,1)*256;
            counter=counter+1;
        end
    end    
    image(:,:,1)=image2;
    image(:,:,2)=image3;
    image(:,:,3)=image4;
    image=uint8(image);
    name=strcat(files{i,1},'-NMF_spikes_22.tif');
    imwrite(image,name,'tif');
end

counter_select=1;
counter=1;
idx=1;
while idx<=numel(sorted_files)
    name=strcat('Maskb_',sorted_files{idx,1});
    A=imread(name);tempidx=unique(A);
    image=zeros(size(A,1),size(A,2),3);
    image2=zeros(size(A));
    image3=zeros(size(A));
    filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\',name(7:end));
    image4=imread(filename);
    image4=single(image4);image4(image4<500)=0;image4=image4/max(max(image4));image4=image4*256;
    j=1;
    while counter <= sorted_index(idx)
        if Index_Select(counter)
            image2(A==tempidx(j))=W_spike_6(counter_select);
            image3(A==tempidx(j))=W_spike_21(counter_select);
            counter=counter+1;
            counter_select=counter_select+1;
        else
            counter=counter+1;
        end
        j=j+1;
    end
    image(:,:,1)=image2;
    image(:,:,2)=image3;
    image(:,:,3)=image4;
    image=uint8(image);
    name=strcat(sorted_files{idx,1},'-NMF_sp_6-21.tif');
    imwrite(image,name,'tif');
    idx=idx+1;
end

counter=1;
counter_select=1;
for i = 1:size(files,1)
	image=zeros(size(Masks{i,1},1),size(Masks{i,1},2),3);
	image2=zeros(size(Masks{i,1}));
    image3=zeros(size(Masks{i,1}));
	name=sorted_files{i,1};
	filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\',name);
    image4=imread(filename);	
    image4=single(image4);image4(image4<500)=0;image4=image4/max(max(image4));image4=image4*256;
    for j = 2:Masks{i,3}
		if Index_Select(counter)
			idx=Masks{i,2}(j);
			image2(Masks{i,1}==idx)=W_spike(counter_select,6);
			image3(Masks{i,1}==idx)=W_spike(counter_select,19);
			counter=counter+1;
			counter_select=counter_select+1;
		else
            counter=counter+1;
        end
    end
	image(:,:,1)=image2;
    image(:,:,2)=image3;
    image(:,:,3)=image4;
    image=uint8(image);
    name=strcat(sorted_files{i,1},'-NMF_sp_6-19.tif');
    imwrite(image,name,'tif');
end
        


counter=1;
counter_select=1;
for i = 1:size(files,1)
	image=zeros(size(Masks{i,1},1),size(Masks{i,1},2),3);
	image2=zeros(size(Masks{i,1}));
    image3=zeros(size(Masks{i,1}));
	image2b=zeros(size(Masks{i,1}));
    image3b=zeros(size(Masks{i,1}));
	name=sorted_files{i,1};
	filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\',name);
    image4=imread(filename);	
    image4=single(image4);image4(image4<500)=0;image4=image4/max(max(image4));image4=image4*256;
    for j = 2:Masks{i,3}
		if Index_Select(counter)
			idx=Masks{i,2}(j);
			image2(Masks{i,1}==idx)=W_ZS(counter_select,26);
			image3(Masks{i,1}==idx)=W_ZS(counter_select,33);
			image2b(Masks{i,1}==idx)=W_ZS(counter_select,26);
			image3b(Masks{i,1}==idx)=W_ZS(counter_select,34);
			counter=counter+1;
			counter_select=counter_select+1;
		else
            counter=counter+1;
        end
    end
	image2=image2/max(max(image2));image2=image2*256;
	image3=image3/max(max(image3));image3=image3*256;
	image2b=image2b/max(max(image2b));image2b=image2b*256;
	image3b=image3b/max(max(image3b));image3b=image3b*256;
	image(:,:,1)=image2;
    image(:,:,2)=image3;
    image(:,:,3)=image4;
    image=uint8(image);
    name=strcat(sorted_files{i,1},'-NMF_ZS_26-33.tif');
    imwrite(image,name,'tif');
	image=zeros(size(Masks{i,1},1),size(Masks{i,1},2),3);
	image(:,:,1)=image2b;
    image(:,:,2)=image3b;
    image(:,:,3)=image4;
    image=uint8(image);
    name=strcat(sorted_files{i,1},'-NMF_ZS_26-34.tif');
    imwrite(image,name,'tif');
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1049, 895]);
imagesc(Hypothesis);colormap(jet);

betas_PCA=zeros(size(aud8freq,1),size(mean_PCA_fluo,1));
for x=1:size(mean_PCA_fluo,1)
    betas_PCA(:,x)=mvregress(aud8freq',mean_PCA_fluo(x,:)');    
end

%aud8freq is the regressor design matrix and mean_NMF are my observations
Hypothesis=zeros(size(betas_PCA,2),size(aud8freq,2));
rsquared=zeros(1,size(betas_PCA,2));
for y=1:size(betas_PCA,2)
    test=zeros(1,size(aud8freq,2));
    for x=1:size(betas_PCA,1)
        temp=betas_PCA(x,y)*aud8freq(x,:);
        test=test+temp;
    end
    Hypothesis(y,:)=test;
    rsquared(y)=sum((Hypothesis(y,:) - mean(mean_PCA_fluo(y,:))).^2) / sum((mean_PCA_fluo(y,:) - mean(mean_PCA_fluo(y,:))).^2);
end

goodbetas=find(rsquared>0.2);
x = linspace(1,size(Hypothesis,2),size(Hypothesis,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1000]);
counter=1;xplot=floor(sqrt(size(goodbetas,2)));yplot=ceil(size(goodbetas,2)/xplot)
for i=goodbetas
    subplot(xplot,yplot,counter);plot(x,mean_PCA_fluo(i,:),x,Hypothesis(i,:));
    counter=counter+1;
end

    
x=1;test=beta(1,x)*aud8freq(1,:);

betas=zeros(size(aud8freq,1),size(Fluo_select,1));
for x=1:size(Fluo_select,1)
betas(:,x)=mvregress(aud8freq',Fluo_select(x,:));
end

counter=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1049, 895]);
x = linspace(1,size(Hypothesis,2),size(Hypothesis,2));
for i=goodbetas
    select{counter}=Fluo_Traces(Index(:,2)==i,:);
    subplot(size(goodbetas,2),1,counter);plot(x,mean(select{counter},1),x,Stimulus);
    counter=counter+1;
end

tic;
distance=zeros(size(DFonF_Dtrend,2),size(DFonF_Dtrend,2));
parfor i=1:size(DFonF_Dtrend,2)
    dist_temp=zeros(1,size(DFonF_Dtrend,2));
    temp=DFonF_Dtrend(:,i);
    for j=1:size(DFonF_Dtrend,2)
        if j>=i
            dist_temp(j)=dtw_c(temp,DFonF_Dtrend(:,j),5);
        end
    end
    distance(i,:)=dist_temp;
end
time=toc;

distance=zeros(size(Fluo,2),size(Fluo,2));
for i=1:size(Fluo,2)
    for j=1:size(Fluo,2)
        distance(i,j)=dtw_c(Fluo(:,i),Fluo(:,j),50);
    end    
end


Kmed_mean=zeros(k,size(Fluo,1));
for i=1:k
    Kmed_mean(i,:) = mean(ZS(:,inds==i),2)';        
end

cluster_mean=zeros(k,size(ZS,2));
for i=1:k
    cluster_mean(i,:) = mean(ZS(T==i,:),1);        
end

R = corrcoef(ZS');
R=triu(R,1);
merged=ZS;
threshold=0.7;
while max(max(R))>threshold
    todelete=[];
    for x=1:size(merged,1)
        tomerge=find(R(x,:)>threshold);
        if ~isempty(tomerge)
            todelete=[todelete tomerge];
            tomerge=[x tomerge];
            merged(x,:)=mean(merged(tomerge,:),1);            
        end
    end
    merged(todelete,:)=[];
    R=corrcoef(merged');
    R=triu(R,1);    
end

NMF_ZS=[];
idx_ZS=[];
mean_NMF_ZS=zeros(size(H_ZS,1),size(H_ZS,2),'double');
for i = 1:size(H_ZS,1)
idx_ZS{i}=find(W_ZS(:,i)>2);
NMF_ZS=ZS(idx_ZS{i},:);
mean_NMF_ZS(i,:)=mean(NMF_ZS,1);
end
clearvars NMF_ZS idx_ZS i;
            
NMF_ZS=[];
idx_ZS=[];
%mean_NMF_sp=zeros(size(H_spike,1),size(H_spike,2),'double');
mean_NMF_fluo=zeros(size(H_ZS,1),size(H_ZS,2),'double');
for i = 1:size(H_ZS,1)
idx_ZS{i}=find(W_ZS(:,i)>5);
%NMF_ZS=Spikes(idx_ZS{i},:);
%mean_NMF_sp(i,:)=mean(NMF_ZS,1);
NMF_ZS=ZS_strong(idx_ZS{i},:);
mean_NMF_fluo(i,:)=mean(NMF_ZS,1);
end
clearvars NMF_ZS idx_ZS i;

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1280, 1024]);
for i =1:3
subplot(3,1,i);scatter3(score(:,1),score(:,2),score(:,3),10,score(:,i),'.');
end

[row col]=find(aud5vol==0.00850739000000000);
corr=zeros(size(Fluo_select,1),1);
%progressbar
parfor x = 1:size(Fluo_select,1)
    for y = 1:numel(col)
        for z = 0:5
            w=col(y)+z;
            temp=Fluo_select(x,w:w+size(GCaMP6,1)-1);
            temp=corrcoef(temp,GCaMP6);
            if temp(1,2)>0.7 %| temp(1,2)<-0.7
                corr(x)=corr(x)+1;
                break
            end            
        end       
    end
    %progressbar(x/size(Fluo_select,1));
end
clearvars x y z w row col;

explained_sum=zeros(size(explained,1),2);
for i = 1:numel(explained)
    explained_sum(i,1)=sum(explained(1:i,:));
    explained_sum(i,2)=sum(explained2(1:i,:));
end

kmean_resp_sp=zeros(max(idx_kmean_sp),size(Responders_sp,2));
for i=1: max(idx_kmean_sp)
    temp=mean(Responders_fluo(idx_kmean_sp==i,:),1);
    kmean_resp_sp(i,:)=temp;
end
clearvars temp;

[row col]=find(aud8freq==0.00850739000000000);
corr=zeros(size(FluoTraces,1),1);
%progressbar
parfor x = 1:size(FluoTraces,1)
    for y = 1:numel(col)
        for z = 0:5
            w=col(y)+z;
            temp=FluoTraces(x,w:w+size(GCaMP6,1)-1);
            temp=corrcoef(temp,GCaMP6);
            if temp(1,2)>0.8 | temp(1,2)<-0.8
                corr(x)=corr(x)+1;
                break
            end            
        end       
    end
    %progressbar(x/size(Fluo_select,1));
end
clearvars x y z w row col;

NMF_ZS=[];
idx_ZS=[];
%mean_NMF_sp=zeros(size(H_spike,1),size(H_spike,2),'double');
%mean_PCA_fluo=zeros(50,size(coeff,2),'double');
mean_PCA_Fluo=zeros(50,size(coeff,2),'double');
for i = 1:size(mean_PCA_Fluo,1)
idx_ZS{i}=find(score_Fluo(:,i)>500);
%NMF_ZS=Spikes(idx_ZS{i},:);
%mean_NMF_sp(i,:)=mean(NMF_ZS,1);
NMF_PCA=Fluo_ZS(idx_ZS{i},:);
mean_PCA_Fluo(i,:)=mean(NMF_PCA,1);
%NMF_PCA=Resp_Fluo(idx_ZS{i},:);
%mean_PCA_fluo(i,:)=mean(NMF_PCA,1);
end
clearvars NMF_ZS NMF_PCA i;

NMF_ZS=[];
idx_HighZS=[];
%mean_NMF_sp=zeros(size(H_spike,1),size(H_spike,2),'double');
mean_PCA_Highfluo=zeros(50,size(coeff,2),'double');
mean_PCA_HighZS=zeros(50,size(coeff,2),'double');
for i = 1:size(mean_PCA_Highfluo,1)
idx_HighZS{i}=find(score2(:,i)>5);
%NMF_ZS=Spikes(idx_ZS{i},:);
%mean_NMF_sp(i,:)=mean(NMF_ZS,1);
NMF_PCA=HighResp_ZS(idx_HighZS{i},:);
mean_PCA_HighZS(i,:)=mean(NMF_PCA,1);
NMF_PCA=HighResp(idx_HighZS{i},:);
mean_PCA_Highfluo(i,:)=mean(NMF_PCA,1);
end
clearvars NMF_ZS NMF_PCA i;

DF=mean(SelectCorr_Dtrend(:,1:50),2);
DF=repmat(DF,1,size(SelectCorr_Dtrend,2));
DFonF_Dtrend=(SelectCorr_Dtrend-DF)./DF;

explained_sum=zeros(size(explained,1),2);
for i = 1:numel(explained)
explained_sum(i,1)=sum(explained(1:i,:));
%explained_sum(i,2)=sum(explained(1:i,:));
end

mean_Kmeans=zeros(max(idxKmeans),size(Sp_infer,2),'double');
for i = 1:max(idxKmeans)
    mean_Kmeans(i,:)=mean(SelectCorr(idxKmeans==i,:),1);
end


mean_Kmeans_Dtrend=zeros(max(idxKmeans_DFDtrend),size(DFonF_Dtrend,2),'double');
for i = 1:max(idxKmeans_DFDtrend)
    mean_Kmeans_Dtrend(i,:)=mean(DFonF_Dtrend(idxKmeans_DFDtrend==i,:),1);
end

betas_Kmeans=zeros(size(aud8freq,1),size(Cmap,1));
for x=1:size(Cmap,1)
    betas_Kmeans(:,x)=mvregress(aud8freq',Cmap(x,:)');    
end

%aud8freq is the regressor design matrix and mean_NMF are my observations
Hypothesis=zeros(size(betas_Kmeans,2),size(aud8freq,2));
rsquared=zeros(1,size(betas_Kmeans,2));
for y=1:size(betas_Kmeans,2)
    test=zeros(1,size(aud8freq,2));
    for x=1:size(betas_Kmeans,1)
        temp=betas_Kmeans(x,y)*aud8freq(x,:);
        test=test+temp;
    end
    Hypothesis(y,:)=test;
    rsquared(y)=sum((Hypothesis(y,:) - mean(mean_Kmeans(y,:))).^2) / sum((mean_Kmeans(y,:) - mean(mean_Kmeans(y,:))).^2);
end

goodbetas=find(rsquared>0.2);
x = linspace(1,size(Hypothesis,2),size(Hypothesis,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1000]);
counter=1;xplot=floor(sqrt(size(goodbetas,2)));yplot=ceil(size(goodbetas,2)/xplot)
for i=goodbetas
    subplot(xplot,yplot,counter);plot(x,mean_Kmeans(i,:),x,Hypothesis(i,:));
    counter=counter+1;
end

for i=1:size(mean_Kmeans,1)
    mdl=stepwiselm(aud8freq',mean_Kmeans(i,:)','linear','Criterion','bic','Upper','interactions');
    model(i).coef=mdl.Coefficients;
    model(i).MSE=mdl.MSE;
    model(i).Fitted=mdl.Fitted;
    model(i).rsquared=mdl.Rsquared.Adjusted;
end

goodbetas=find([model.rsquared]>0.3);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1000]);
counter=1;xplot=floor(sqrt(size(goodbetas,2)));yplot=ceil(size(goodbetas,2)/xplot)
for i=goodbetas
    subplot(xplot,yplot,counter);plot(x,mean_Kmeans(i,:),x,model(i).Fitted);
    counter=counter+1;
end

opol = 6;
mean_Kmeans_detrend=zeros(size(mean_Kmeans));
t = (1:size(mean_Kmeans,2))';
for i=1:size(mean_Kmeans,1)
    [p,s,mu] = polyfit(t,mean_Kmeans(i,:)',opol);
    f_y = polyval(p,t,[],mu);
    mean_Kmeans_detrend(i,:)=mean_Kmeans(i,:)-f_y';
end

for i=1:size(Cmap_BF2,2)
    mdl=stepwiselm(aud8freq',Cmap_BF3(:,i),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions');
    model3_baseline(i).coef=mdl.Coefficients;
    model3_baseline(i).MSE=mdl.MSE;
    model3_baseline(i).Fitted=mdl.Fitted;
    model3_baseline(i).rsquared=mdl.Rsquared.Adjusted;
end

goodbetas_baseline=find([model3_baseline.rsquared]>0.5);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1000]);
counter=1;xplot=floor(sqrt(size(goodbetas_baseline,2)));yplot=ceil(size(goodbetas_baseline,2)/xplot);
for i=goodbetas_baseline
    subplot(xplot,yplot,counter);plot(x,Cmap_BF3(:,i),x,model3_baseline(i).Fitted);
    counter=counter+1;
end


opol = 5;
SelectCorr_Dtrend=zeros(size(SelectCorr));
t = (1:size(SelectCorr,2))';
for i=1:size(SelectCorr,1)
    [p,s,mu] = polyfit(t,SelectCorr(i,:)',opol);
    f_y = polyval(p,t,[],mu);
    p=mean(SelectCorr(i,1:50),2);
    SelectCorr_Dtrend(i,:)=SelectCorr(i,:)-f_y';
    SelectCorr_Dtrend(i,:)=SelectCorr_Dtrend(i,:)+p;
end
clearvars opol t f_y s p mu

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1000]);
subplot(2,1,1);imagesc(mean_Kmeans,[0 0.3]);colormap jet
subplot(2,1,2);imagesc(mean_Kmeans_Dtrend,[0 0.2]);colormap jet

Index_temp1=find(Index_Select==1);
Index_temp2=find(Nb_corr>3);
Index_final=Index_temp1(Index_temp2);

pts=[1 18 210 420 670];
[mean_Kmeans_Dtrend_BF,yfit] = bf(mean_Kmeans_Dtrend',pts,'linear');

temp=[];
for i=1:length(goodbetas_Dtrend)
    test=find(idxKmeans_Dtrend==i);
    temp{i}=Index_final(test);
end
Index_good=zeros(sum(cellfun('length',temp)),2);
Index_good(1:length(temp{1}),1)=temp{1};Index_good(1:length(temp{1}),2)=1;
counter=length(temp{1});
for i=2:length(goodbetas_Dtrend)
    Index_good(counter+1:counter+length(temp{i}),1)=temp{i};Index_good(counter+1:counter+length(temp{i}),2)=i;
    counter=counter+length(temp{i});
end
clearvars test temp

index_nogang=zeros(length(sorted_files),1);
for i = 1:size(sorted_files,1)    
    name=strcat('Maskb_',sorted_files{i,1});    
    A=imread(name);
    index_nogang(i)=max(max(A));    
end
clearvars A i

index_nogang_sum=zeros(length(index_nogang),1);
for i = 1:length(index_nogang_sum)
    index_nogang_sum(i)=sum(index_nogang(1:i));
end
clearvars i

test=zeros(length(SelectCorr_nogang),1);
match=zeros(length(SelectCorr),2);
for i = 889:size(SelectCorr,1)
    test=zeros(length(SelectCorr_nogang),1);
    for j = 1:size(SelectCorr_nogang,1)
        temp=corrcoef(SelectCorr(i,1:150),SelectCorr_nogang(j,1:150));
        test(j)=temp(1,2);        
        if temp(1,2)>0.98
            break
    end
    [M,I]=max(test);
    match(i,1)=M;match(i,2)=I;
end
clearvars test temp M I