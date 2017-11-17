
x = inputdlg('Enter Number of repeats', 'Number of repeats');
NbRepeats = str2num(x{:});

listOfNMFmat = dir('NMF-*.mat');
numberOfNMF = numel(listOfNMFmat);
mydata(numberOfNMF).H=0;
mydata(numberOfNMF).W=0;
for j = 1:numberOfNMF
    try
    myfilename = sprintf('NMF-%d.mat',j);
    mydata(j)=load(myfilename);
    end
end


for j = 1:NbRepeats:numberOfNMF
    F=[];
    G=[];
    %try
   for k=0:(NbRepeats-1)
    F=horzcat(F,mydata(j+k).H);
    G=vertcat(G,mydata(j+k).W);
   end
    [cr,lgs] = xcorr(F,0,'coeff');
    test2 = find(cr>0.8);
    start=0;
    temp=0;
    cr2=reshape(cr,NbRepeats*size(mydata(1).H,2),NbRepeats*size(mydata(1).H,2));
    allcr{j}=cr2;    
    f=figure('Visible', 'off');imagesc(triu(cr2)>0.8); axis image;
    myfilename = sprintf('cross-correlation-%d.png',j);
    saveas(f, myfilename,'png');
    clear temparray;
    clear H;
    clear I;
    for idx = 1:numel(test2)
        element = test2(idx);
        x = ceil(element/(NbRepeats*size(mydata(1).H,2)));
        y = mod(element,NbRepeats*size(mydata(1).H,2));
        if y == 0
            y=NbRepeats*size(mydata(1).H,2);
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
    H=vertcat(H,temparray);
    Hval{j}=I;
    Wval{j}=H;
    %end
end

save('NMF - collapse.mat','Hval','Wval','-v7.3')
save('cross-correlation - NMF.mat','allcr','-v7.3')
close all;

for i=1:size(Hval,2)
	if isempty(Hval{i})
	else
		for k=1:size(Hval{i},2)
			f=figure('Visible', 'off');subplot(1,2,1);plot(Hval{i}(:,k));subplot(1,2,2);imagesc(squeeze(Wval{i}(k,:,:)),[0,3]); axis image;axis off;
			myfilename = sprintf('NMFcollapse-%d-%d.png',i,k);
			saveas(f, myfilename,'png');
            close all;
		end
	end
end

%Tone after vis
stim=zeros(5,500);
stim(1,85:100)=1;
stim(2,130:145)=1;
stim(3,175:190)=1;
stim(4,220:235)=1;
stim(5,305:315)=1;
stim(5,325:335)=1;
stim(5,350:360)=1;
stim(5,370:380)=1;
stim(5,395:415)=1;

%Tone with vis
stim2=zeros(5,500);
stim2(1,85:100)=1;
stim2(2,130:145)=1;
stim2(3,175:190)=1;
stim2(4,220:235)=1;
stim2(5,85:95)=1;

%Tone before vis
stim4=zeros(5,500);
stim4(1,285:300)=1;
stim4(2,330:345)=1;
stim4(3,375:390)=1;
stim4(4,420:435)=1;
stim4(5,55:65)=1;
stim4(5,75:85)=1;
stim4(5,100:110)=1;
stim4(5,120:130)=1;
stim4(5,145:155)=1;

%WJ
stim5=zeros(6,500);
stim5(1,185:200)=1;
stim5(2,230:245)=1;
stim5(3,275:290)=1;
stim5(4,320:335)=1;
stim5(5,25:35)=1;
stim5(5,45:55)=1;
stim5(5,70:80)=1;
stim5(5,90:100)=1;
stim5(5,115:125)=1;
stim5(6,405:415)=1;
stim5(6,455:465)=1;

%WJ2
stim6=zeros(7,500);
stim6(1,185:205)=1;
stim6(2,230:250)=1;
stim6(3,275:295)=1;
stim6(4,320:340)=1;
stim6(5,25:40)=1;
stim6(5,45:60)=1;
stim6(5,70:85)=1;
stim6(5,90:105)=1;
stim6(5,115:130)=1;
stim6(6,405:425)=1;
stim6(7,455:475)=1;

%5tones
stim7=zeros(7,600);
stim7(1,20:35)=1;
stim7(1,200:215)=1;
stim7(1,400:415)=1;
stim7(1,20:35)=1;
stim7(2,45:60)=1;
stim7(2,225:240)=1;
stim7(2,425:440)=1;
stim7(3,70:85)=1;
stim7(3,250:265)=1;
stim7(3,450:465)=1;
stim7(4,95:110)=1;
stim7(4,275:290)=1;
stim7(4,475:490)=1;
stim7(5,120:135)=1;
stim7(5,300:315)=1;
stim7(5,500:515)=1;
stim7(6,145:160)=1;
stim7(6,325:340)=1;
stim7(6,525:540)=1;
stim7(7,170:185)=1;
stim7(7,350:365)=1;
stim7(7,550:565)=1;

%VisSBIb-c
stim8=zeros(8,400);
stim8(1,55:70)=1;
stim8(2,100:115)=1;
stim8(3,145:160)=1;
stim8(4,190:205)=1;
stim8(5,205:220)=1;
stim8(6,250:265)=1;
stim8(7,295:310)=1;
stim8(8,340:355)=1;

%Tone_vis_WJ
stim9=zeros(10,600);
stim9(1,55:65)=1;
stim9(1,73:83)=1;
stim9(1,98:108)=1;
stim9(2,160:170)=1;
stim9(3,200:225)=1;
stim9(4,250:275)=1;
stim9(5,290:315)=1;
stim9(6,303:323)=1;
stim9(7,350:370)=1;
stim9(8,395:420)=1;
stim9(9,440:465)=1;
stim9(10,505:530)=1;
stim9(10,555:580)=1;

%Tone_vis_WJ_grav
toneviswjgrav=zeros(6,600);
toneviswjgrav(1,52:61)=1;
toneviswjgrav(1,72:80)=1;
toneviswjgrav(1,96:105)=1;
toneviswjgrav(2,165:180)=1;
toneviswjgrav(3,210:225)=1;
toneviswjgrav(4,250:265)=1;
toneviswjgrav(5,352:367)=1;
toneviswjgrav(5,402:417)=1;
toneviswjgrav(6,452:467)=1;
toneviswjgrav(6,502:517)=1;

%SLM test Kevin
SLMKev=zeros(1,200);
SLMKev(1,25:35)=1;
SLMKev(1,88:98)=1;
SLMKev(1,148:158)=1;

%aud8freq
    aud8freq=zeros(8,700);
    aud8freq(1,21:21+5)=1;
    aud8freq(1,406:406+5)=1;
    aud8freq(1,611:611+5)=1;

    aud8freq(2,46:46+5)=1;
    aud8freq(2,381:381+5)=1;
    aud8freq(2,536:536+5)=1;

    aud8freq(3,71:71+5)=1;
    aud8freq(3,356:356+5)=1;
    aud8freq(3,511:511+5)=1;

    aud8freq(4,96:96+5)=1;
    aud8freq(4,331:251);
    aud8freq(4,561:561+5)=1;

    aud8freq(5,121:121+5)=1;
    aud8freq(5,306:306+5)=1;
    aud8freq(5,486:486+5)=1;

    aud8freq(6,146:146+5)=1;
    aud8freq(6,281:281+5)=1;
    aud8freq(6,586:586+5)=1;

    aud8freq(7,171:171+5)=1;
    aud8freq(7,256:256+5)=1;
    aud8freq(7,461:461+5)=1;

    aud8freq(8,196:196+5)=1;
    aud8freq(8,231:231+5)=1;
    aud8freq(8,636:636+5)=1;

stim3=aud8freq;

for i=1:size(Hval,2)
	if isempty(Hval{i})
    else
        RHO=[];
		ccf=[];
		pval=[];
		for k=1:size(Hval{i},2)
			for j=1:size(stim3,1)
				[RHO,P]=corrcoef(Hval{i}(:,k),transpose(stim3(j,:)));
				ccf(j,k)=RHO(1,2);
				pval(j,k)=P(1,2);
			end
		end
		f=figure('Visible', 'off');imagesc(ccf);axis image;
		myfilename = sprintf('NMFcorrelation-%d.png',i);
		saveas(f, myfilename,'png');
        close all;
		test2=find(ccf>0.5);
		for idx = 1:numel(test2)
			element = test2(idx);
			y = ceil(element/size(stim3,1));
			x = mod(element,size(stim3,1));
            if x==0
                x=size(stim3,1);
            end
			f=figure('Visible', 'off');subplot(1,2,1);plot(Hval{i}(:,y));subplot(1,2,2);imagesc(squeeze(Wval{i}(y,:,:)),[0,3]); axis image;axis off;
			myfilename = sprintf('NMF-%d_%d_stim-%d.png',i,y,x);
			saveas(f, myfilename,'png');
            close all;
		end
	end
end

