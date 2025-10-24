%Read RNS data

[dex,dat]=ReadECoGData('C:\Users\u244026\Documents\HaneefRNS\132783014100930000.dat', 'test.lay');
dat=dat{1};
dat=dat';
[xdn,~,~,spike]=dbtDenoise(dat, 250, 0.1);


[b,a]=butter(4,[4 8]/125,'bandpass');
filter(b,a,xdn);
th=filter(b,a,xdn);
th_power=bandpower(th);

[b,a]=butter(4,[8 12]/125,'bandpass');
filter(b,a,xdn);
alph=filter(b,a,xdn);
alph_power=bandpower(alph);

[b,a]=butter(4,[12 30]/125,'bandpass');
filter(b,a,xdn);
bet=filter(b,a,xdn);
bet_power=bandpower(bet);

[b,a]=butter(4,[30 70]/125,'bandpass');
filter(b,a,xdn);
gam=filter(b,a,xdn);
gam_power=bandpower(gam);

[b,a]=butter(4,[70 124.5]/125,'bandpass');
filter(b,a,xdn);
hgam=filter(b,a,xdn);
hgam_power=bandpower(hgam);

figure,bar([th_power,alph_power,bet_power,gam_power,hgam_power]);

phas=angle(hilbert(th));
amp=  abs(hilbert(gam));
amph= abs(hilbert(hgam));

[MIlo,KL]=modulationIndex(phas,amp)
[MIhi,KL2]=modulationIndex(phas,amph)