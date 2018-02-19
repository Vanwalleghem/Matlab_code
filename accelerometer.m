mean=[];
std=[];
for i=0:7
    idx=8-i;      
    mean(idx)=nanmean(norm_filt(86800-(i*10250):88200-(i*10250)));
    std(idx)=nanstd(norm_filt(86800-(i*10250):88200-(i*10250)));
end
g=9.80665;
mean=mean*g;
std=std*g;

Fstop = 70;
Fpass = 85;
Astop = 80;
Apass = 1;
Fs = 2e3;

% d = designfilt('highpassfir','StopbandFrequency',Fstop, ...
%   'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
%   'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','equiripple');
% fvtool(d)

d = designfilt('highpassiir','StopbandFrequency',Fstop ,...
  'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
  'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','butter');

fvtool(d2)

dataIn = eight_freq(:,2:4); dataOut = filter(d,dataIn);figure;plot(dataOut(:,1))

Fstop = 300;
Fpass = 350;
Astop = 65;
Apass = 0.5;
Fs = 2e3;

% d = designfilt('highpassfir','StopbandFrequency',Fstop, ...
%   'PassbandFrequency',Fpass,'StopbandAttenuation',Astop, ...
%   'PassbandRipple',Apass,'SampleRate',Fs,'DesignMethod','equiripple');
% fvtool(d)

d2 = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',340,'HalfPowerFrequency2',440, ...
    'SampleRate',2000);

volumeOut = filter(d2,five_volumes(:,2:3));%figure;plot(volumeOut(:,1));hold on;plot(five_volumes(:,2));
volumeOut = volumeOut*9.81;

d3 = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',70,'HalfPowerFrequency2',120, ...
    'SampleRate',2000);

dataOut_1 = filter(d3,dataIn);figure;plot(dataOut_1(:,1))

[Xa Ya] = alignsignals(dataOut(:,1), dataOut(:,2));
Temp=[Xa(1:length(Ya)) Ya];
[Xa Ya] = alignsignals(Temp(:,1), dataOut(:,3));
find (Xa==0);figure;plot(Temp);hold on;plot(Ya(2000:length(Temp)));
Temp2=Ya(2000:length(Ya));Temp2=Temp2*9.81;
Temp=Temp*9.81;

freq_acc=[];
start=5600;freq_acc(1,1)=max(Temp(start:start+1000,1));freq_acc(1,2)=max(Temp(start:start+1000,2));freq_acc(1,3)=max(Temp2(start:start+1000));
i=2;start=15500;freq_acc(i,1)=max(Temp(start:start+1000,1));freq_acc(i,2)=max(Temp(start:start+1000,2));freq_acc(i,3)=max(Temp2(start:start+1000));
i=3;start=26000;freq_acc(i,1)=max(Temp(start:start+1000,1));freq_acc(i,2)=max(Temp(start:start+1000,2));freq_acc(i,3)=max(Temp2(start:start+1000));
i=4;start=36000;freq_acc(i,1)=max(Temp(start:start+1000,1));freq_acc(i,2)=max(Temp(start:start+1000,2));freq_acc(i,3)=max(Temp2(start:start+1000));
i=5;start=46500;freq_acc(i,1)=max(Temp(start:start+1000,1));freq_acc(i,2)=max(Temp(start:start+1000,2));freq_acc(i,3)=max(Temp2(start:start+1000));
i=6;start=56700;freq_acc(i,1)=max(Temp(start:start+1000,1));freq_acc(i,2)=max(Temp(start:start+1000,2));freq_acc(i,3)=max(Temp2(start:start+1000));
i=7;start=66700;freq_acc(i,1)=max(Temp(start:start+1000,1));freq_acc(i,2)=max(Temp(start:start+1000,2));freq_acc(i,3)=max(Temp2(start:start+1000));
i=8;start=77200;freq_acc(i,1)=max(Temp(start:start+1000,1));freq_acc(i,2)=max(Temp(start:start+1000,2));freq_acc(i,3)=max(Temp2(start:start+1000));

[Xa Ya] = alignsignals(volumeOut(:,1), volumeOut(:,2));
volumeOut=[Xa Ya(1:length(Xa))];

vol_acc=[];
start=22500;vol_acc(1,1)=max(volumeOut(start:start+1000,1));vol_acc(1,2)=max(volumeOut(start:start+1000,2));
i=2;start=32800;vol_acc(i,1)=max(volumeOut(start:start+1000,1));vol_acc(i,2)=max(volumeOut(start:start+1000,2));
i=3;start=43040;vol_acc(i,1)=max(volumeOut(start:start+1000,1));vol_acc(i,2)=max(volumeOut(start:start+1000,2));
i=4;start=53200;vol_acc(i,1)=max(volumeOut(start:start+1000,1));vol_acc(i,2)=max(volumeOut(start:start+1000,2));
i=5;start=63200;vol_acc(i,1)=max(volumeOut(start:start+1000,1));vol_acc(i,2)=max(volumeOut(start:start+1000,2));
vol_accb=vol_acc*2;