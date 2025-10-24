ptnum='_7532636';
ptInit='AEH';

pathP='C:\Users\u244026\Documents\HaneefRNS\Hippo\';
fold=['Baylor College of Medicine, Houston_',ptInit,ptnum,' EXTERNAL #PHI']
subfold=['Baylor College of Medicine, Houston_' ptInit,ptnum, ' Data EXTERNAL #PHI']


Sheet=importdata(fullfile(pathP, fold,'\', 'Baylor College of Medicine, Houston_AEH_7532636_ECoG_Catalog.csv')) 
SheetP=Sheet.textdata;

%add on extra column (23) to give original index of event 
for i=1:length(SheetP)
    Sheetmat=1:length(SheetP);
    SheetP{i,23}={Sheetmat(i)};
end

Sheetdex=find(strcmp(SheetP(:,8),'Scheduled')==1);
SheetS=SheetP(Sheetdex,:);

daydex=zeros(length(SheetS),1);
for i=1:length(SheetS)
    timeS=SheetS{i,4};
    timeD=str2num(timeS(12:13));
   if timeD <= 20 && timeD >=08 ==1
       daydex(i)=1;
   end
end
daydex=logical(daydex);
SheetD=SheetS(daydex,:);

%emptycount=0;
fullcount=0;
th_power=nan(length(SheetD),1);
alph_power = nan(length(SheetD), 1);
bet_power = nan(length(SheetD), 1);
gam_power = nan(length(SheetD), 1);
hgam_power = nan(length(SheetD), 1);
MIlo=nan(length(SheetD),1);
MIhi=nan(length(SheetD),1);
MIthresh=nan(length(SheetD),1);
hMIthresh=nan(length(SheetD),1);


for i=1:length(SheetD)
datfile=SheetD{i,6};

datlocp=fullfile(pathP, fold, '\', subfold);
foldDir=dir(datlocp);
Searchdex=find(contains({foldDir.name},datfile)==1);

if isempty(Searchdex)
 %   emptycount=emptycount+1;

else
   fullcount=fullcount+1;
%i
%nono=1
datloc=fullfile(pathP, fold, '\', subfold, '\', datfile);
[dex,dat]=ReadECoGData(datloc,  'test.lay');

for q=1:length(dat)
q

datchan=dat{q};
datchan=datchan';

[xdn,~,~,spike]=dbtDenoise(datchan, 250, 0.1);


[b,a]=butter(4,[4 8]/125,'bandpass');
filter(b,a,xdn);
th=filter(b,a,xdn);
th_power(i)=bandpower(th);

[b,a]=butter(4,[8 12]/125,'bandpass');
filter(b,a,xdn);
alph=filter(b,a,xdn);
alph_power(i)=bandpower(alph);

[b,a]=butter(4,[12 30]/125,'bandpass');
filter(b,a,xdn);
bet=filter(b,a,xdn);
bet_power(i)=bandpower(bet);

[b,a]=butter(4,[30 70]/125,'bandpass');
filter(b,a,xdn);
gam=filter(b,a,xdn);
gam_power(i)=bandpower(gam);

[b,a]=butter(4,[70 124.5]/125,'bandpass');
filter(b,a,xdn);
hgam=filter(b,a,xdn);
hgam_power(i)=bandpower(hgam);

%figure,bar([th_power,alph_power,bet_power,gam_power,hgam_power]);

phas=angle(hilbert(th));
amp=  abs(hilbert(gam));
amph= abs(hilbert(hgam));

[MIlo(i),KL]=modulationIndex(phas,amp);
[MIhi(i),KL2]=modulationIndex(phas,amph);


%MI shuffle to create normalized distribution
perm=250;
MI_R=zeros(perm,1);
hMI_R=zeros(perm,1);

for p=1:perm
shift=randi([0 length(phas)]);
phasR=circshift(phas,shift);

MI_R(p)=modulationIndex(phasR, amp);
[MImean, sigma]=normfit(MI_R);
MIthresh(i)=norminv(0.95, MImean,sigma);

hMI_R(p)=modulationIndex(phasR, amph);
[hMImean, hsigma]=normfit(hMI_R);
hMIthresh(i)=norminv(0.95, hMImean,hsigma);
end


end
end

nandex=~isnan(th_power); %filter out absent data files
th_powerR=th_power(~isnan(th_power));
alph_powerR=alph_power(~isnan(alph_power));
bet_powerR = bet_power(~isnan(bet_power));
gam_powerR = gam_power(~isnan(gam_power));
hgam_powerR = hgam_power(~isnan(hgam_power));
MIlo_Obs = MIlo(~isnan(MIlo));  %observed MI
MIhi_Obs = MIhi(~isnan(MIhi));
MIthresh=MIthresh(nandex);  %95 percentile MI of chance distribution
hMIthresh=hMIthresh(nandex);

numPAClo=sum(MIloR>MIthresh); %number of PAC instances
numPAChi=sum(MIhiR>hMIthresh);

%save values of interest into spreadsheet

HippodatRNS{1,1}='channelnum';
HippodatRNS{1,2}='channelID';%.channelID=strcat(string(q),ptnum);
HippodatRNS{1,3}='channelname'; %=SheetP{2,11+q}
% Save the results into a structured array for further analysis
HippodatRNS{1,4}='th_power'; %= th_powerR;
HippodatRNS{1,5}='alph_power'; % = alph_powerR;
HippodatRNS{1,6}='bet_power'; 
HippodatRNS{1,7}='gam_power';% = gam_powerR;
HippodatRNS{1,8}='hgam_power';% = hgam_powerR;
HippodatRNS{1,9}='MIlo_Obs';% = MIlo_Obs;
HippodatRNS{1,10}='MIhi_Obs';% = MIhi_Obs;
HippodatRNS{1,11}='numPAClo'; % = numPAClo;
HippodatRNS{1,12}='numPAChi'; % = numPAChi;
HippodatRNS{1,13}='MIthresh'; %=MIthresh;
HippodatRNS{1,14}='hMIthresh'; %=MIthresh;

HippodatRNS{q+1,1}=q;
HippodatRNS{q+1,2}=strcat(string(q),ptnum);
HippodatRNS{q+1,3}=SheetP{2,11+q}
% Save the results into a structured array for further analysis
HippodatRNS{q+1,4}= th_powerR;
HippodatRNS{q+1,5}= alph_powerR;
HippodatRNS{q+1,6} = bet_powerR;
HippodatRNS{q+1,7} = gam_powerR;
HippodatRNS{q+1,8} = hgam_powerR;
HippodatRNS{q+1,9} = MIlo_Obs;
HippodatRNS{q+1,10} = MIhi_Obs;
HippodatRNS{q+1,11} = numPAClo;
HippodatRNS{q+1,12} = numPAChi;
HippodatRNS{q+1,13}= MIthresh;
HippodatRNS{q+1, 14}=hMIthresh;

end
% Save the structured array  for future analysis
save('HippodatRNS_results_dummy.mat', 'HippodatRNS');









