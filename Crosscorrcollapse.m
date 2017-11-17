x = inputdlg('Enter Number of repeats', 'Number of repeats');
NbRepeats = str2num(x{:});
listOfcorrmat = dir('correlation-*-0.mat');
listOfcorr = dir('correlation-1-*.mat');
numberofcorr = numel(listOfcorr);
numberOfcorrmat = numel(listOfcorrmat);
%mycorr(numberOfcorrmat).correlation0=0;

for j = 1:numberOfcorrmat
	for i = 0:(numberofcorr-1)
		myfilename = sprintf('correlation-%d-%d.mat',j,i);
		myvariable = sprintf('correlation%d',i);
		temp=load(myfilename,myvariable);
		mycorr{j,i+1}=temp.(myvariable);
    end
end
Gcorr=[];

for i = 1:numberofcorr
	for j = 1:NbRepeats:numberOfcorrmat
        Fcorr=cat(3,mycorr{j:j+(NbRepeats-1),i});
		Gcorr{ceil(j/NbRepeats),i}=mean(Fcorr,3);
	end
end

save('avgCorr.mat','Gcorr','-v7.3')

for i=1:size(Gcorr,1)
	for j=1:size(Gcorr,2)
		if isempty(Gcorr{i,j})
		else
			f=figure('Visible', 'off'); imagesc(Gcorr{i,j},[-0.3,0.7]); axis image;axis off;
			myfilename = sprintf('correlation-%d-%d.png',i,j-1);
			saveTightFigure(f, myfilename);
			myfilename = sprintf('correlation-%d-%d.fig',i,j-1);
			saveTightFigure(f, myfilename);
			close all;
		end
	end
end
