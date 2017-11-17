clear all;    

[fName,fDir,fIndex] = uigetfile(...
{'*.tif','TIFF Image Stack (*.tif)';'*.mat','Matlab Data File (*.mat)'},...
pwd);
fFile = fullfile(fDir,fName);

ds = 1;
A = imfinfo(fFile);
I = cell(length(A),1);

for i = 1:length(A)
Inw = imread(fFile,i);
if (size(Inw,3) == 3)
else
I{i} = downsample(downsample(double(Inw),ds)',ds)';
end

if (i == 1)
[Imx,Imn,Isum] = deal(I{i});
else
[Imx,Imn,Isum] = deal(max(Imx,I{i}),...
min(Imn,I{i}),Isum+I{i});
end
end


Imean = Isum/length(I);
im = Imean;
imMean = mean(mean(im));
im2 = im;
clims(1) = prctile(reshape(im,1,[]),0.15);
clims(2) = prctile(reshape(im,1,[]),99.5);
J = adapthisteq(uint16(im),'NumTiles',[ceil(size(im,1)/64) ceil(size(im,2)/64)],'ClipLimit', 0.02,'NBins', 1024,'Distribution','uniform');figure;imagesc(J);axis image;
level = graythresh(J)/2;
im16=uint16(J);
imTempMask = im2bw(im16,level);
mask = medfilt2(imTempMask,[6,6]);


% Filtering Image prior to Segmentation



%figure; imagesc(J); colormap gray; axis image;

Fxy=double(J);
smwindow1=4;
smwindow2=2;
smwindow3=6;
SE = strel('disk',3,0);
filter1=fspecial('log', [smwindow1 smwindow1], smwindow1/4);
filter2=fspecial('gaussian', [smwindow2 smwindow2], smwindow2/4);
filter3=fspecial('gaussian', [smwindow3 smwindow3], smwindow3/4);

SmFresc1 = imfilter(Fxy,filter2);
SmFxy=imfilter(SmFresc1,filter1);
SmFxyInv = SmFxy*(-1);
SmFxyInvFilt = imfilter(SmFxy,filter3);
SmFxyTHInv = imtophat(SmFxyInvFilt, SE);
SmFresc1 = imfilter(SmFxyTHInv,filter2);

% %             % Segmentation
L = watershed(SmFresc1.*mask);

STATS = regionprops(L, 'Area','PixelIdxList','Eccentricity','Centroid','Extrema','Solidity');

wneuron = zeros(size(STATS,1),1);
for i = 1:size(STATS,1)
    if STATS(i).Area(1,1)>0
        wneuron(i,1)=1;
    end
end
roilist = find(wneuron);


%             Collect fluorescence data for each ROI over time
data = zeros(size(I,1),size(roilist,1));
for j = 1:size(I,1)
    for i = 1:size(roilist,1)
        data(j,i) = nanmean(I{j,1}(STATS(roilist(i)).PixelIdxList));
    end
end

lpFilt = designfilt('lowpassfir', 'FilterOrder', 2,...
    'PassbandFrequency', 1, 'StopbandFrequency', 2, ...
    'SampleRate', 5, 'DesignMethod', 'ls');
data2 = filtfilt(lpFilt,data);


% Calculate Correlated ROIs
m = size(data,1);
data3 = data2-repmat(mean(data2),size(data2,1),1);
data4sd = std(data3);    
Rho = zeros(size(data,2),size(data,2));
for i = 1:size(roilist,1)
    Sigma = (data3(:,i)' * data3) / (m-1);
    Rho(i,:) = Sigma ./ (data4sd(1,i)'*data4sd);
end
highRho = find(Rho>=0.99);
IND = [];
[IND(:,1),IND(:,2)] = ind2sub(size(Rho),highRho);
abcd = IND(:,1)-IND(:,2); abcde = find(abcd>0);
correlatedROIs = IND(abcde,:);
correlatedROIs2 = []; correlatedROIs3 = [];
correlatedROIs2(:,1) = roilist(correlatedROIs(:,1));
correlatedROIs2(:,2) = roilist(correlatedROIs(:,2));
correlatedROIs3 = sortrows(correlatedROIs2,1);
correlatedROIs4 = correlatedROIs3;
ss = var(data2);
for i = 1:size(correlatedROIs3,1)
    if ss(1,correlatedROIs3(i,1))>=500
        correlatedROIs4(i,:) = NaN;
    end
end
correlatedROIs5 = correlatedROIs4(isnan(correlatedROIs4(:,1))==0,:);
correlatedROIs6 = unique(correlatedROIs5(:,1));


% Merge adjacent correlated ROIs
SE2 = strel('square',2);
SE3 = strel('square',2);
Ltemp = logical(L);
loop = 1;

for j = correlatedROIs6(:,1)'
    temtemp = find(correlatedROIs5==j);
    temtemp = temtemp(temtemp<size(correlatedROIs5,1));
    temtemp2 = correlatedROIs5(temtemp,2);
    temp1 = ismember(L,[j;temtemp2]);
    temp2 = imdilate(temp1,SE2);
    temp3 = imerode(temp2,SE3);
    temp4 = []; temp4(:,:,1) = temp3; temp4(:,:,2) = Ltemp; 
    temp5 = max(temp4,[],3);
    Ltemp = temp5;
    loop = loop+1;
end    

XL = watershed((Ltemp*-1));
STATS_temp = regionprops(XL, 'Area','PixelIdxList','Eccentricity','Centroid','Extrema','Solidity');

areamax = 180;
areamin = 20;
excentmax = 0.9;
solidmin = 0.75;

a=field2num(STATS_temp,'Area');
ex=field2num(STATS_temp,'Eccentricity');
sol = field2num(STATS_temp,'Solidity');
wneuron1=ones(size(a,1),1);
w=find(a>areamax);
wneuron1(w)=0;
w=find(a<areamin);
wneuron1(w)=0;
w=find(ex>excentmax);
wneuron1(w)=0;
w=find(sol<solidmin);
wneuron1(w)=0;
Nneurons=squeeze(sum(wneuron1));
neuronlist = find(wneuron1);

[Lia, Locb] = ismember(XL, neuronlist);
Locb2 = Locb;
dLia = double(Lia);
figure; imagesc(dLia); colormap(gray); axis image;

%USE THIS TO GET A LOOK AT THE SEGMENTED ROIs OVER RAW IMAGE

im2 = im;
gtmax = gt(im2,clims(2));
ltmin = lt(im2,clims(1));
im_sc = im2;
im_sc(gtmax) = clims(2);
im_sc(ltmin) = clims(1);
im3 = double(im2).*dLia;
im4 = double(im_sc)-(im3*0.08);
clear rgb;
rgb(:,:,1) = im4*1.8;
rgb(:,:,2) = im4*1.8;
rgb(:,:,3) = im2*1;
rgb = uint16(rgb);
%imwrite(rgb,[fName 'segment.tif'],'tif');
RGB64 = double(rgb)/65535;
rgblims(1) = prctile(reshape(RGB64,1,[]),0.2);
rgblims(2) = prctile(reshape(RGB64,1,[]),99.5);
JJJ = imadjust(RGB64,rgblims);
figure; imagesc(JJJ); axis image;




% Collect fluorescence data for each ROI over time
data = zeros(size(I,1),size(neuronlist,1));
for j = 1:size(I,1)
    for i = 1:size(neuronlist,1)
        data(j,i) = nanmean(I{j,1}(STATS_temp(neuronlist(i)).PixelIdxList));
    end
end


Fbackground = prctile(reshape(im,1,[]),60);
Fbackground = repmat(Fbackground,numel(I),size(data,2));
data2 = data - Fbackground;
Fnaught = prctile(data2([1:9,31:40],:),25);
Fnaught = repmat(Fnaught,numel(I),1);
deltaFonF = 100.*((data2-Fnaught)./Fnaught);





