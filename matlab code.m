

F=horzcat(A.H, B.H, C.H, D.H, E.H)
[cr,lgs] = xcorr(F,5,'coeff');
test2 = find(cr(6,:)>0.8)
test = cr > 0.8
G=vertcat(A.W, B.W, C.W, D.W, E.W)
test4=reshape(test3,45,45)
cr2=reshape(squeeze(cr(6,:)),45,45)
cr2=triu(cr2)


start=0
temp=0
clear temparray;
clear H;
clear I;
for idx = 1:numel(test2)
	element = test2(idx);
	x = ceil(element/27);
	y = mod(element,27);
	if y == 0
		y=27;
	end
	if any(temp == y)
	else
		temp(end+1) = y;
		if start==x
			if exist('temparray') == 1
				temparray=vertcat(temparray,G(y,:,:));
			else
				temparray=G(y,:,:);
			end
		else
			if exist('H') == 1
				temparray = mean(temparray,1);
				H=vertcat(H,temparray);
				temparray=G(x,:,:);
			else
				if exist('temparray') == 1
					temparray = mean(temparray,1);
					H=temparray;
					temparray=G(x,:,:);
				else
					temparray=G(x,:,:);
				end
			end
			if exist('I') == 1
				I=horzcat(I,F(:,x));
			else
				I=F(:,x);
			end
		end
	end
	start=x;
end


figure;
for i=1:size(H,1)
	if i==size(H,1)
		subplot(5,4,i); imagesc(squeeze(H(1,:,:)),[0,3]); axis image;axis off;
		sub_pos = get(gca,'position');
	else
		subplot(5,4,i); imagesc(squeeze(H(i+1,:,:)),[0,3]); axis image;axis off; 
		sub_pos = get(gca,'position'); % get subplot axis position
	end
set(gca,'position',sub_pos.*[0.6 1 1.2 1.2]) % stretch its width and height
colormap jet;
end

figure;
for i=1:size(I,2)
subplot(5,4,i); plot(I(:,i));
end


for i=1:size(H,1)
	if i==size(H,1)
		figure; imagesc(squeeze(H(1,:,:)),[0,3]); axis image;axis off;
	else
		figure; imagesc(squeeze(H(i+1,:,:)),[0,3]); axis image;axis off;
	end
end
for i=1:size(I,2)
figure; plot(I(:,i));
end



j=1;
figure; 
for i=1:size(H,1)
	if i==size(H,1)
		subplot(7,4,j); imagesc(squeeze(H(1,:,:)),[0,3]); axis image;axis off;
		subplot(7,4,j+1);plot(I(:,i)); axis on;
		sub_pos = get(gca,'position');
	else
		subplot(7,4,j); imagesc(squeeze(H(i+1,:,:)),[0,3]); axis image;axis off; 
		subplot(7,4,j+1);plot(I(:,i)); axis on;
		sub_pos = get(gca,'position'); % get subplot axis position
	end
	%set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
	colormap jet;
	j=j+2;
end

>> stim(2,130:145)=1;
>> stim(3,175:190)=1;
>> stim(4,220:235)=1;
>> stim(5,55:65)=1;
>> stim(5,75:85)=1;
>> stim(5,100:110)=1;
>> stim(5,120:130)=1;
>> stim(5,145:155)=1;

RHO=[]
for i=1:size(I,2)
	for j=1:size(stim,2)
		[RHO,P]=corrcoef(I(:,i),stim(:,j));
		if exist('ccf') == 1
				ccf=vertcat(ccf,RHO);
				pval=vertcat(pval,P);
			else
				ccf=RHO;
				pval=P;
			end
	end
end

RHO=[]
ccf=[]
pval=[]
for i=1:size(I,2)
	for j=1:size(stim,2)
		[RHO,P]=corrcoef(I(:,i),stim(:,j));
		ccf(i,j)=RHO(1,2);
		pval(i,j)=P(1,2);
	end
end


for i=1:size(I,2)
	[RHO,P]=corrcoef(I(:,i),stim(:,1));
	if exist('ccf') == 1
				ccf=vertcat(ccf,RHO);
				pval=vertcat(pval,P);
			else
				ccf=RHO;
				pval=P;
			end
	end



115 - 194

for j= 13:18 
	for i= 0:5
		for k = 1:5
			if (j==39 && k==3) || (j==39 && k==2 && i==5)
			elseif (j==23 && k==5)
			elseif j<=36
			myfilename = sprintf('correlation-%d-%d.mat',(j-1)*5+k,i);
			myvarname = sprintf('correlation%d',i);
			mydata(j).(myvarname)(k) = load(myfilename,myvarname);
			else
			myfilename = sprintf('correlation-%d-%d.mat',(j-1)*5+k+1,i);
			myvarname = sprintf('correlation%d',i);
			mydata(j).(myvarname)(k) = load(myfilename,myvarname);
			end
		end
	end
end

F=[];
avg=struct;
for j= 13:18 
	for i = 0:5
		myvarname = sprintf('correlation%d',i);
		for k = 1:5
			if (j==39 && k==3) || (j==39 && k==2 && i==5)
			elseif (j==23 && k==5)
			else
			F=cat(3,F,mydata(j).(myvarname)(k).(myvarname));
			end
		end
	G=mean(F,3);
	avg(j).(myvarname)=G;
	F=[];
	end
end

for j= 13:18
	figure;
	for i=1:6
		myvarname = sprintf('correlation%d',i-1);
		subplot(2,3,i);imagesc(avg(j).(myvarname),[0,0.8]);axis image;axis off;
		sub_pos = get(gca,'position'); % get subplot axis position
		set(gca,'position',sub_pos.*[1 1 1.2 1.2]) % stretch its width and height
		colormap jet;
	end
end

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
G=[];

for i = 1:numberofcorr
	for j = 1:NbRepeats:numberOfcorrmat
		F=cat(3,mycorr{j,i},mycorr{j+1,i},mycorr{j+2,i});
		G{ceil(j/NbRepeats),i}=mean(F,3);
	end
end

for i=1:size(G,1)
	for j=1:size(G,2)
		if isempty(G{i,j})
		else
			f=figure('Visible', 'off'); imagesc(G{i,j},[-0.3,0.7]); axis image;axis off;
			myfilename = sprintf('correlation-%d-%d.png',i,j-1);
			saveTightFigure(f, myfilename);
			myfilename = sprintf('correlation-%d-%d.fig',i,j-1);
			saveTightFigure(f, myfilename);
			close all;
		end
	end
end
