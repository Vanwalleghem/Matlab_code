listOflinreg = dir('*f-linreg.mat');
numberOflinreg = numel(listOflinreg);
x=load(listOflinreg(1).name);
numberofbetas=size(x.betas,2);

for j = 1:numberOflinreg
	x=load(listOflinreg(j).name);
    name=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\Masks\Mask_',listOflinreg(j).name(1:end-11));
    A=imread(name);
    image=zeros(size(A,1)*size(A,2),3);
    for i = 1:size(x.betas,1)
        if x.rsq(i)>0.2
            if (size(find(A==i),1)<120) && (8<size(find(A==i),1))
                image(A==i,1)=x.betas(i,1);
                image(A==i,2)=x.betas(i,2);
                image(A==i,3)=x.betas(i,3);
            end
        end
    end
    image=reshape(image,size(A,1),size(A,2),3);
    name=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\',listOflinreg(j).name(1:end-11));
    B=imread(name);
    C=imfuse(image,B,'blend');
    name=strcat('img_spikes_',listOflinreg(j).name(1:end-14),'png');
    imwrite(C,name);
end

image=zeros(size(A,1)*size(A,2),3);
for i = 1:size(betas,1)
image(A==i,1)=betas(i,1);
image(A==i,2)=betas(i,2);
image(A==i,3)=betas(i,3);
end
image=reshape(image,size(A,1),size(A,2),3);

imshow( image2 ); hold on;
h = imagesc( Bnorm ); % show the edge image
set( h, 'AlphaData', .8 ); % .5 transparency
colormap gray

counter=1;
for i = 1:size(files,1)
    name=strcat('Mask_',files_sort{i,1});
    ispresent = cellfun(@(s) ~isempty(strfind(name,s)),listOflinregb);
    ispresent=find(ispresent(1,:)==1);
    image=zeros(size(A{ispresent,1}));
    for j = 1:max(max(A{ispresent,1}))
        image(A{ispresent,1}==j)=W_hat(counter,22);
        counter=counter+1;
    end
    image=uint16(image);
    name=strcat(files_sort{i,1},'-NMF_spikes_22.tif');
    imwrite(image,name,'tif');
end

listOfMasks = dir('Maskb_*.tif');
listOfMasks=struct2cell(listOfMasks);
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


W_ZS_33=W_ZS_33/norm(W_ZS_33,Inf);
counter=1;
for i=1:size(Filelistb,2)
    image=zeros(size(A3{i,1},1),size(A3{i,1},2),3);
    image2=zeros(size(A3{i,1}));    
    filename=strcat('D:\Pictures\processed\Tonotropy\GCaMP6f - 8 tones\',Filelistb{1,i});
    image4=imread(filename);
    image4=single(image4);image4(image4<500)=0;image4=image4/max(max(image4));image4=image4*256;
    for j = 1:max(max(size(A3{i,1})))
        if counter < index(i,1)
            image2(A3{i,1}==j)=W_ZS_33(counter,1)*256;
            counter=counter+1;
        end
    end
    image(:,:,1)=image2;
    image(:,:,3)=image4;
    image=uint8(image);
    name=strcat(files{i,1},'-NMF_spikes_33.tif');
    imwrite(image,name,'tif');
end
