x = inputdlg('Enter Number of repeats', 'Number of repeats');
NbRepeats = str2num(x{:});
listOflinreg = dir('linreg-*.mat');
numberOflinreg = numel(listOflinreg);
load('linreg-1.mat');
numberofbetas=size(betas,1);

for j = 1:numberOflinreg
	myfilename = sprintf('linreg-%d.mat',j);
    x=load(myfilename,'betas');
	mylinreg{j}=x.betas;	
end
Glinr=[];

for j = 1:NbRepeats:numberOflinreg
    Flinr=cat(4,mylinreg{1,j:j+(NbRepeats-1)});
	Glinr{ceil(j/NbRepeats)}=mean(Flinr,4);
end

save('avgLinR.mat','Glinr','-v7.3')

for j=1:size(Glinr,2)
    f=figure('Visible', 'off');
    for i=1:numberofbetas
        if isempty(Glinr{1,j})
        else
            subplot(floor(sqrt(numberofbetas)),ceil(numberofbetas/floor(sqrt(numberofbetas))),i); imagesc(squeeze(Glinr{1,j}(i,:,:))); axis image;axis off;
        end
    end
    f=tightfig(f);
   myfilename = sprintf('linreg_avg-%d.png',j);
   saveas(f, myfilename);
   myfilename = sprintf('linreg_avg-%d.fig',j);
   saveas(f, myfilename);
   close all;
end
