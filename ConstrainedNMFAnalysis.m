MatFiles=dir('*.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
Spikes=load(name, 'Spikes');
Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);
for i = 2:length(MatFiles)
name=strcat(MatFiles(i).name);
C=load(name, 'DenoisedTraces');
C=C.DenoisedTraces;
    if i==3
        C=[C(:,1) C(:,1) C(:,1:58)];
    end
S=load(name, 'Spikes');
S=S.Spikes;
N=load(name, 'Noise');
N=N.Noise;
F=load(name, 'idx_components');
F=F.idx_components+1;
GC=C(F,:);
Noise=vertcat(Noise,N);
Calcium=vertcat(Calcium,C);
Spikes=vertcat(Spikes,S);
Fitness=horzcat(Fitness,F);
GoodCalcium=vertcat(GoodCalcium,GC);
MatFiles(i).number=size(Calcium,1);
end
clearvars GC C S F N name i;

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
