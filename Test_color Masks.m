[sorted_index,idxsort_files]=sortrows(index);
sorted_files=files(idxsort_files);
A=[];
for i=1:length(sorted_files)
    filename=strcat('Maskb_',sorted_files{i});
    A{i}=imread(filename);    
end
clearvars i idxsort_files;

sorted_index2=zeros(size(sorted_index));
for i=1:length(A)
    sorted_index2(i)=max(max(A{i}));
end
clearvars i filename

sorted_index3=zeros(size(sorted_index));
sorted_index3(1)=double(max(max(A{1})));
for i=2:length(A)    
    sorted_index3(i)=sorted_index3(i-1)+double(max(max(A{i})));
end
clearvars i

counter=1;
idx_corr=Nb_corr>3;
idx_final=zeros(length(idx_corr),1);
for i=1:length(idx_corr)
    if idx_corr(i)
        idx_final(i)=idxKmeans_DF(counter);
        counter=counter+1;
    end
end
clearvars counter;

counter=1;
rsquared_final=zeros(length(idx_corr),1);
for i=1:length(idx_corr)
    if idx_corr(i)
        rsquared_final(i)=rsquared(counter);
        counter=counter+1;
    end
end
clearvars counter;


colors = distinguishable_colors(length(SelectSClusters_BF_05));
colors = colors*256;
for i=1:length(A)
    B=A{i};B=double(B);B(B==0)=NaN;
    image=zeros(size(B,1),size(B,2),3);
    for j=min(min(B)):max(max(B))
        if any(idx_final(j)==SelectSClusters_BF_05)
            clust=find(idx_final(j)==SelectSClusters_BF_05);
            for k=1:3              
              image2=zeros(size(B));
              image2(B==j)=colors(clust,k);
              image(:,:,k)=image(:,:,k)+image2;
            end
        end
    end
    image=uint8(image); 
    name=strcat('05rsq_',sorted_files{i},'-Kmeans.tif');
    imwrite(image,name,'tif');
end
clearvars image image2 i j k name B clust;

files_F8=dir('05rsq_*F8*Kmeans.tif');
for i = 1:length(files_F8)
    name=strcat(files_F8(i).name);
    B=imread(name);         
    filename=strcat('D:\Pictures\processed\Tonotropy\Tectum\',name(7:length(name)-11));
    image4=imread(filename);image4=single(image4);image4(image4<500)=0;image4=image4/max(max(image4));image4=image4*256;    
    image=double(B)+(double(repmat(image4,1,1,3)));
    image=uint8(image); 
    name=strcat('BackGD_Kmeans05_',name(7:length(name)-11));
    imwrite(image,name,'tif');        
end
clearvars image image4 i name B files_F8 filename;

counter=1;m=length(Model_all);
progressbar;
colors = [[127,201,127],[190,174,212],[253,192,134],[255,255,153],[56,108,176],[240,2,127],[191,91,23],[102,102,102]];
colors = reshape(colors,3,8);colors = colors';
for i=1:length(A)
    B=A{i};B=double(B);B(B==0)=NaN;maxrsq=0;
    image=zeros(size(B,1),size(B,2),3);        
    for j=min(min(B)):max(max(B))
        temp=zeros(size(B));
        if idx_corr(j)
            for k=1:8            
                if Model_all(counter).coef.pValue(k) < 0.01
                    for rgb=1:3
                        temp(B==j)=Model_all(counter).rsquared*colors(k,rgb)/8;
                        image(:,:,rgb)=image(:,:,rgb)+temp;
                        if maxrsq<Model_all(counter).rsquared
                            maxrsq=Model_all(counter).rsquared;
                        end
                    end
                end
            end
        counter=counter+1;
        end
        progressbar(counter/m)
    end
    image=image/max(max(max(image)));image=image*maxrsq*256;
    image=uint8(image); 
    name=strcat('all_rsqColor_',sorted_files{i});
    imwrite(image,name,'tif');    
end
clearvars image i j name B counter m colors maxrsq;

for i=1:length(A)
    B=A{i};B=double(B);B(B==0)=NaN;
    image=zeros(size(B,1),size(B,2),3);        
    for j=min(min(B)):max(max(B))
        temp=zeros(size(B));      
        for k=4:6
            if Model_all(j).coef.pValue(k) < 0.01              
                temp(B==j)=Model_all(j).coef.Estimate(k)*Model_all(j).rsquared;
                image(:,:,k-3)=image(:,:,k-3)+temp;
            end
        end
    end
    image=image/max(max(max(image)));image=image*256;
    image=uint8(image); 
    name=strcat('rsqColor2_',sorted_files{i});
    imwrite(image,name,'tif');
end
clearvars image i j name B;

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


files_F8=dir('rsqColor_*F8*.tif');
for i = 1:length(files_F8)
    name=strcat(files_F8(i).name);
    B=imread(name);         
    filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\AVG_',name(10:length(name)));
    image4=imread(filename);image4=single(image4);image4(image4<500)=0;image4=image4/max(max(image4));image4=image4*256;    
    image=double(B)+(double(repmat(image4,1,1,3)));
    image=uint8(image); 
    name=strcat('BackGD_rsqColor_',name(10:length(name)));
    imwrite(image,name,'tif');        
end
clearvars image image4 i name B files_F8 filename;

rsquared=[Model_all.rsquared]';
counter=1;
rsquared_final=zeros(length(idx_corr),1);
for i=1:length(idx_corr)
    if idx_corr(i)
        rsquared_final(i)=rsquared(counter);
        counter=counter+1;
    end
end
clearvars counter rsquared;

for i=1:length(A)
B=A{i};B=double(B);B(B==0)=NaN;
image=zeros(size(B));
for j=min(min(B)):max(max(B))
image(B==j)=rsquared_final(j)*256;
end
image=uint8(image);
name=strcat('rsq_',sorted_files{i});
imwrite(image,name,'tif');
end
clearvars image i j name B;

for i=1:length(A)
    B=A{i};B=double(B);B(B==0)=NaN;
    image=zeros(size(B,1),size(B,2),3);
    for j=min(min(B)):max(max(B))
        if any(idxKmeans_ZS(j)==SelectClusters_ZS)
            clust=find(idxKmeans_ZS(j)==SelectClusters_ZS);
            for k=1:3              
              image2=zeros(size(B));
              image2(B==j)=colors(clust,k)*256;
              image(:,:,k)=image(:,:,k)+image2;
            end
        end
    end
    image=uint8(image); 
    name=strcat('Kmeansmodel',sorted_files{i},'-Kmeans.tif');
    imwrite(image,name,'tif');
end
clearvars image image2 i j k name B clust;

counter=1;m=length(model);
progressbar;
colors = [[256,0,0],[0,0,256]];
colors = reshape(colors,3,2);colors = colors';
for i=1:length(A)
    B=A{i};B=double(B);B(B==0)=NaN;maxrsq=0;
    image=zeros(size(B,1),size(B,2),3);        
    for j=min(min(B)):max(max(B))
        temp=zeros(size(B));
        if length(model(counter).coef.pValue)>1
            for k=1:1+(length(model(counter).coef.pValue)-2)
                if model(counter).coef.pValue(k+1) < 0.01
                    for rgb=1:3
                        temp(B==j)=model(counter).rsquared*colors(k,rgb);
                        image(:,:,rgb)=image(:,:,rgb)+temp;
                        if maxrsq<model(counter).rsquared
                            maxrsq=model(counter).rsquared;
                        end
                    end
                end
            end
        end
        counter=counter+1;        
        progressbar(counter/m)
    end
    image=uint8(image); 
    name=strcat('all_rsqColor_',sorted_files{i});
    imwrite(image,name,'tif');    
end
clearvars image i j name B counter m colors maxrsq;

%colors = distinguishable_colors(length(Goodbetas_BF));
colors = [[1 0 0]; [0 1 0]; [0 0 1]];
colors = colors*256;
Cluster_map=[];
for i=1:length(A)
    B=A{i};B=double(B);B(B==0)=NaN;
    filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\AVG_',sorted_files{i});
    image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(max((image4)));image4=image4*200;image4=double(repmat(image4,1,1,3));image4=uint8(image4);
    %image=zeros(size(B,1),size(B,2),3);
    for j=min(min(B)):max(max(B))
        if any(idx_final(j)==Goodbetas_Select)
            clust=find(idx_final(j)==Goodbetas_Select);
            for k=1:3              
              image2=image4(:,:,k);
              image2(B==j)=colors(clust,k);
              image4(:,:,k)=image2;
            end
        end
    end
    Cluster_map{i}=
    %image=uint8(image); 
    name=strcat('Final_',sorted_files{i});
    %filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\AVG_',sorted_files{i});
    %image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(prctile(image4,80));image4=image4*200;
    %image4=double(repmat(image4,1,1,3));image4=uint8(image4)+image;
    image4=uint8(image4);
    imwrite(image4,name,'tif');
end
clearvars image image2 image4 i j k name filename B clust;

for i = 1:length(files_F8)
    name=strcat(files_F8(i).name);
    B=imread(name);         
    filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\AVG_',name(10:length(name)));
    image4=imread(filename);image4=single(image4);image4(image4<500)=0;image4=image4/max(prctile(image4,80));image4=image4*200;    
    image=double(B)+(double(repmat(image4,1,1,3)));
    image=uint8(image); 
    name=strcat('BackGD_rsqColor_',name(10:length(name)));
    imwrite(image,name,'tif');        
end
clearvars image image4 i name B files_F8 filename;

colors = colors*256;
for i=1:length(A)
    B=A{i};B=double(B);B(B==0)=NaN;
    image=zeros(size(B,1),size(B,2),3);
    for j=min(min(B)):max(max(B))
        if idx_final(j)>0            
            for k=1:3              
              image2=zeros(size(B));
              image2(B==j)=colors(idx_final(j),k);
              image(:,:,k)=image(:,:,k)+image2;
            end
        end
    end
    image=uint8(image); 
    name=strcat('Model_Kmeans_',sorted_files{i});
    filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\AVG_',sorted_files{i});
    image4=imread(filename);image4=double(image4);image4(image4<500)=0;image4=image4/max(prctile(image4,80));image4=image4*200;
    image4=double(repmat(image4,1,1,3));image4=uint8(image4)+image;
    imwrite(image4,name,'tif');
end
clearvars image image2 image4 i j k name filename B clust;

Index_Select
counter=1;counter2=1;
idx_corr=Nb_corr>3;
idx_final=zeros(length(Index_Select),1);
for i=1:length(Index_Select)
    if Index_Select(i)        
        if idx_corr(counter)
            idx_final(i)=idxKmeans_DF_20(counter2);
            counter2=counter2+1;
        counter=counter+1;
    end
end
clearvars counter;

counter=1;counter2=1;
idx_final=zeros(length(idx_corr),1);
for i=1:length(idx_corr)
    if idx_corr(i)
        if any(counter==idx_thr5)
            counter=counter+1;
        else
            idx_final(i)=idxKmeans_DF_20(counter2);
            counter=counter+1;counter2=counter2+1;
        end
    end
end
clearvars counter counter2;
