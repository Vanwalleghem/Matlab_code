x = inputdlg('Enter Number of repeats', 'Number of repeats');
NbRepeats = str2num(x{:});
listOfcorrmat = dir('correlation-*.mat');
numberOfcorrmat = numel(listOfcorrmat);
temp=load('correlation-1.mat');
numberofcorr=size(temp.correlations,1);

for j = 1:numberOfcorrmat
	myfilename = sprintf('correlation-%d.mat',j);
    x=load(myfilename,'correlations');
	mycorr{j}=x.correlations;	
end
Gcorr=[];

for j = 1:NbRepeats:numberOfcorrmat
    Fcorr=cat(4,mycorr{1,j:j+(NbRepeats-1)});
	Gcorr{ceil(j/NbRepeats)}=mean(Fcorr,4);
end

save('avgCorr.mat','Gcorr','-v7.3')

for j=1:size(Gcorr,2)
    f=figure('Visible', 'off');
    for i=1:numberofcorr
        if isempty(Gcorr{1,j})
        else
            subplot(floor(sqrt(numberofcorr)),ceil(numberofcorr/floor(sqrt(numberofcorr))),i); imagesc(squeeze(Gcorr{1,j}(i,:,:))); axis image;axis off;
        end
    end
    f=tightfig(f);
   myfilename = sprintf('correlation_avg-%d.png',j);
   saveas(f, myfilename);
   close all;
end

