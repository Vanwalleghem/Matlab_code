listOfmat = dir('*-Mask representation.mat');
numberOfmat = numel(listOfmat);
x=load(listOfmat(1).name);
goodNMFs=[6,7,9,11,16,17,21,27];

for j = 1:numberOfmat    
    x=load(listOfmat(j).name);    
    h = figure;set(h, 'Visible', 'off');
    k=1;
    for ii = goodNMFs;
        subaxis(4,2,k, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0); imagesc(squeeze(x.NMF(ii,:,:)));set(gca, 'XTick', []);set(gca, 'YTick', []);axis tight; axis off;axis image;k=k+1; 
        %subplot(4,2,k); imagesc(squeeze(x.NMF(ii,:,:)));set(gca, 'XTick', []);set(gca, 'YTick', []); axis image;k=k+1; 
    end    
    myfilename = regexp(listOfmat(j).name,'([G]\w+-\w+-[F]\d_\d+um)','match');
    saveas(h, myfilename{1},'png');    
end

for j = 1:numberOfmat    
    x=load(listOfmat(j).name);
    myfilename = regexp(listOfmat(j).name,'([G]\w+-\w+-[F]\d_\d+um)','match');
    filepath=strcat('D:\Pictures\processed\20150904\',myfilename{1},'.tif');
    y=imread(filepath,1);
    h = figure;set(h, 'Visible', 'off');
    a = normalization(x.NMF(6,:,:));a(2,:,:)=normalization(x.NMF(7,:,:));a(3,:,:)=normalization(x.NMF(9,:,:));
    a=permute(a,[2,3,1]);
    b = normalization(x.NMF(11,:,:));b(2,:,:)=normalization(x.NMF(16,:,:));b(3,:,:)=normalization(x.NMF(17,:,:));
    b=permute(b,[2,3,1]);
    c = normalization(x.NMF(21,:,:));c(2,:,:)=normalization(x.NMF(27,:,:));c(3,:,:)=zeros(size(x.NMF(9,:,:)));
    c=permute(c,[2,3,1]);
    subaxis(2,3,1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);d = bar(H(:,[6:7,9]));d(1).FaceColor='r';d(2).FaceColor='g';d(3).FaceColor='b';
    d(1).EdgeColor='r';d(2).EdgeColor='g';d(3).EdgeColor='b';
    subaxis(2,3,2, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);e = bar(H(:,[11,16,17]));e(1).FaceColor='r';e(2).FaceColor='g';e(3).FaceColor='b';
    e(1).EdgeColor='r';e(2).EdgeColor='g';e(3).EdgeColor='b';
    subaxis(2,3,3, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);f = bar(H(:,[21,27]));f(1).FaceColor='r';f(2).FaceColor='g';
    f(1).EdgeColor='r';f(2).EdgeColor='g';
    %subaxis(3,2,k, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0); imagesc(squeeze(x.NMF(ii,:,:)));set(gca, 'XTick', []);set(gca, 'YTick', []);axis tight; axis off;axis image; 
    %subplot(4,2,k); imagesc(squeeze(x.NMF(ii,:,:)));set(gca, 'XTick', []);set(gca, 'YTick', []); axis image;k=k+1; 
    subaxis(2,3,4, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);imshow(a);hold on;h=imshow(y,[300 6000]);set(h,'AlphaData',0.3);hold off;
    subaxis(2,3,5, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);imshow(b);hold on;h=imshow(y,[300 6000]);set(h,'AlphaData',0.3);hold off;
    subaxis(2,3,6, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);imshow(c);hold on;h=imshow(y,[300 6000]);set(h,'AlphaData',0.3);hold off;
    
    saveas(h, myfilename{1},'png');
    close all;
end


for j = 1:numberOfmat    
    x=load(listOfmat(j).name);
    myfilename = regexp(listOfmat(j).name,'([G]\w+-\w+-[F]\d_\d+um)','match');
    filepath=strcat('D:\Pictures\processed\20150904\',myfilename{1},'.tif');
    y=imread(filepath,1);
    h = figure;set(h, 'Visible', 'off');
    a = normalization(x.betas(1,:,:));a(2,:,:)=normalization(x.betas(2,:,:));a(3,:,:)=normalization(x.betas(3,:,:));
    a=permute(a,[2,3,1]);
    b = normalization(x.betas(4,:,:));b(2,:,:)=normalization(x.betas(5,:,:));b(3,:,:)=normalization(x.betas(6,:,:));
    b=permute(b,[2,3,1]);
    c = normalization(x.betas(7,:,:));c(2,:,:)=normalization(x.betas(8,:,:));c(3,:,:)=zeros(size(x.betas(1,:,:)));
    c=permute(c,[2,3,1]);
    subaxis(1,3,1, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);imshow(a);hold on;h=imshow(y,[300 6000]);set(h,'AlphaData',0.3);hold off;
    subaxis(1,3,2, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);imshow(b);hold on;h=imshow(y,[300 6000]);set(h,'AlphaData',0.3);hold off;
    subaxis(1,3,3, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);imshow(c);hold on;h=imshow(y,[300 6000]);set(h,'AlphaData',0.3);hold off;
    filename=strcat(myfilename{1},'-betas');
    saveas(h, filename,'png');
    close all;
end
