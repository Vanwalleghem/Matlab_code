filelist=dir('*.tif');
for File=1:length(filelist)
    [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(filelist(File).name,[],50,1,'C:\Temp\CROP');
    [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, [1:size(mixedsig,1)],  0.8,size(mixedsig,1),randn(size(mixedsig,1), size(mixedsig,1)),1e-6,50000);
    PCA_ICA_results(File).ica_sig=ica_sig;
    PCA_ICA_results(File).ica_filters=ica_filters;
end
save('PCA_ICA_results.mat');

figure;
File=28;
Y=bigread2(filelist(File).name,1,100);
Y=mean(Y,3);
CellsortICAplot('contour', PCA_ICA_results(File).ica_filters, PCA_ICA_results(File).ica_sig, Y, [], 0.2, 3, 1);

CellsortICAplot(mode, ica_filters, ica_sig, f0, tlims, dt, ratebin, plottype, ICuse, spt, spc)
CellsortICAplot('contour', ica_filters, ica_sig, f0, [0 300], 0.2, 3, 1);
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, [], 0.8);
[ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, [1:20], 0.8);
CellsortICAplot('contour', ica_filters, ica_sig, f1_d1, [0 300], 0.2, 3, 1);
CellsortICAplot('contour', ica_filters, ica_sig, f1_d1, [], 0.2, 3, 1);
[ica_segments, segmentlabel, segcentroid] = CellsortSegmentation(ica_filters, smwidth, thresh, arealims, plotting)
[ica_segments, segmentlabel, segcentroid] = CellsortSegmentation(ica_filters, 6, 3, 8, 0);
CellsortICAplot('contour', ica_filters, ica_sig, f1_d1, [], 0.2, 3, 1);
load('D:\Downloads\stim.mat')
figure;imagesc(STIM);
figure;imagesc([STIM STIM]);
CellsortICAplot('contour', ica_filters, ica_sig, f1_d1, [0 1500], 0.2, 3, 1);
CellsortICAplot('contour', ica_filters, ica_sig, f1_d1, [0 300], 0.2, 3, 1);

flow=zeros(6,655);
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
%GCaMP6=[0.000256990000000000;0.00850739000000000;0.0654158300000000;0.0784609000000000;0.0764130100000000;0.0665958600000000;0.0579028900000000;0.0467942900000000;0.0232079800000000;0.0144564400000000;0.00695772000000000;0.00526551000000000;0.00299500000000000;0.00198520000000000;0.00128512000000000;0.00134175000000000;0.000403170000000000;0];
back=[57 257 457];
back_off=[106 306 506];
fwd=[157 357 557];
fwd_off=[207 407 607];
flow(1,back(1):back(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(1,back(2):back(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(1,back(3):back(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(1):back_off(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(2):back_off(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(3):back_off(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(1):fwd(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(2):fwd(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(3):fwd(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(1):fwd_off(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(2):fwd_off(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(3):fwd_off(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(5,back(1):back(1)+43)=1;
flow(5,back(2):back(2)+43)=1;
flow(5,back(3):back(3)+43)=1;
flow(6,fwd(1):fwd(1)+43)=1;
flow(6,fwd(2):fwd(2)+43)=1;
flow(6,fwd(3):fwd(3)+43)=1;
clearvars GCaMP6 back back_off fwd fwd_off;


for idx=1:length(PCA_ICA_results)
    ModelResults=[];
    for trace=1:size(PCA_ICA_results(idx).ica_sig,1)
        mdl=stepwiselm(flow',PCA_ICA_results(idx).ica_sig(trace,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
        ModelResults(trace)=mdl.Rsquared.Adjusted;
    end
    PCA_ICA_results(idx).rsquared=ModelResults;
end

for idx=17:length(PCA_ICA_results)
    ModelResults=[];
    for trace=1:size(PCA_ICA_results(idx).ica_sig,1)
        mdl=stepwiselm(flow',PCA_ICA_results(idx).ica_sig(trace,50:704),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
        ModelResults(trace)=mdl.Rsquared.Adjusted;
    end
    PCA_ICA_results(idx).rsquared=ModelResults;
end

idx=1;AllTraces=PCA_ICA_results(idx).ica_sig;
for idx=1:length(PCA_ICA_results)
%     if idx>16
%         AllTraces=vertcat(AllTraces,PCA_ICA_results(idx).ica_sig(:,50:704));
%     else
        AllTraces=vertcat(AllTraces,PCA_ICA_results(idx).ica_sig);
    end
end

for idx=1:length(PCA_ICA_results)
    good_rsq=find(PCA_ICA_results(idx).rsquared>0.2);
    PCA_ICA_results(idx).goodTraces=PCA_ICA_results(idx).ica_sig(good_rsq,:);
end

GoodTraces=PCA_ICA_results(1).goodTraces;
for idx=1:length(PCA_ICA_results)
    if PCA_ICA_results(idx).goodTraces
        if idx>16
            GoodTraces=vertcat(GoodTraces,PCA_ICA_results(idx).goodTraces(:,50:704));
        else
            GoodTraces=vertcat(GoodTraces,PCA_ICA_results(idx).goodTraces);
        end
    end
end

options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(GoodTraces,5,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');

GoodBetas_select=[1 2 3 4 5];
x = linspace(1,size(Cmap_ZS,2),size(Cmap_ZS,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=GoodBetas_select    
    NumberOfCells=length(find(idxKmeans_ZS==i));
    %subplot(5,1,counter);plot(x,Cmap_ZS(i,:),x,Model_DF(i).Fitted);title(num2str(NumberOfCells))
    subplot(xplot,yplot,counter);plot(Cmap_ZS(i,:));title(num2str(NumberOfCells))
    %subplot(xplot,yplot,counter);imagesc(DF(find(idxKmeans_ZS==i),:),[0 0.3]);colormap hot;title(num2str(NumberOfCells))
    xlim([0 size(Cmap_ZS,2)])
    counter=counter+1;
end

ZS=zscore(AllTraces,1,2);
options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(AllTraces,10,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
